#!/usr/bin/env bash

# =====                         =====
# ===== define helper functions ===== 
# =====                         =====

# usage info to stdout and exit 1
usage() {
    echo ""
    echo "Usage with reads:"
    echo "  $0 -1 <read1 fastq> -2 <read2 fastq> -d <refdir> [-r <run label>]"
    echo "Usage with simulation:"
    echo "  $0 -s <sim refs tsv> -d <refdir> [-n <num pairs>] [-r <run label>]"
    echo ""
    echo "  -s  TSV (or CSV) with refpath : percentage (int) newline delimited pairs for simulation. Simulated fastqs will {percent}% of reads from {refpath}"
    echo "  -1  Path to R1 FASTQ file. Both -1 and -2 must be used, the pair holds precedence over -s"
    echo "  -2  Path to R2 FASTQ file. Both -1 and -2 must be used, the pair holds precedence over -s"
    echo "  -d  Directory containing reference directories (each containing the ref FASTA file)"
    echo "  -n  Number of pairs to simulate (positive integer, only valid with -s)"
    echo "  -r  Label for current run, appended to output dir pathname. Default pathname is timestamp (YYYY-MM-DD_hh-mm-ss)"
    echo "  -h  Show this help message"
    echo ""
    exit 1
}

# get line number of file, excludes last empty newline if present, for slurm array batching
count_lines() {
    local filepath=$1
    local count_lines=$(wc -l < "$filepath")

    if [[ $(tail -c1 "$filepath") == $'\n' ]]; then
        count_lines=$(( count_lines - 1 ))
    fi

    echo $count_lines
}

# wait for slurm job to complete before pipeline can continue
await_job_completion() {
    local jobid=$1
    while squeue | col1 | grep $jobid > /dev/null 2>&1; do
        sleep 1
    done
}

[[ $# -eq 0 ]] && usage  # display usage on no flag run


# =====                     =====
# ===== Parse CLI arguments =====
# =====                     =====


while getopts ":1:2:d:s:n:r:h" opt; do
    case $opt in
        1) read1="$OPTARG" ;;
        2) read2="$OPTARG" ;;
        d) refdir="${OPTARG%%/}" ;;
        s) simrefs="$OPTARG" ;;
        n) numpairs="$OPTARG" ;;
        r) runname="$OPTARG" ;;
        h) usage ;;
        \?) echo "Error: Invalid option -$OPTARG" >&2; usage ;;
        :) echo "Error: Option -$OPTARG requires an argument." >&2; usage ;;
    esac
done

# =====                                  =====
# ===== More initialization + validation =====
# =====                                  =====

FLAG_SIMULATION_MODE=0

# Require -d
if [[ ! -d "$refdir" ]]; then
    echo "Error: -d <refdir> must be a valid directory." >&2
    usage
fi

if [[ ! -f "$read1" || ! -f "$read2" ]]; then

    if [[ -f "$simrefs" ]]; then
        echo "Queued: Simulating reads..."
        FLAG_SIMULATION_MODE=1
    else
        echo "Error: Must specify either -1/-2 pair OR -s <sim refs tsv>." >&2
        usage
    fi

    if [[ -z "$numpairs" ]]; then
        numpairs=1000000
    elif ! [[ "$numpairs" =~ ^[0-9]+$ ]] || (( numpairs <= 0 )); then
        echo "Error: -n must be a positive integer." >&2
        usage
    fi
else
    echo "Status: Using passed reads"
fi

outdir="results/$(date +%F_%H-%M-%S)"

if [[ -n "$runname" ]]; then
    outdir="${outdir}-${runname}"
    echo "Status: output directory: $outdir"
fi


PROJ="/gpfs/commons/home/apandit/gghome/lions"
TEMP="/nfs/scratch/$outdir"
TEMPSHARE="/nfs/scratch/results/shared"
BIN="$PROJ/bin"

mkdir -p $outdir
mkdir -p $TEMP
mkdir -p $TEMPSHARE

# =====                     =====
# ===== Actual Script Logic =====
# =====                     =====

#
# 1. download univec
#

echo "Status: Downloading univec..."
curl -Ls -o $TEMP/UniVec_Core.fa ftp://ftp.ncbi.nlm.nih.gov/pub/UniVec/UniVec_Core

#
# 2. cache large files in temp dir
#

non_indexed_refs=$TEMP/non-indexed-references.txt
touch $non_indexed_refs

echo "Status: Copying uncached references to temp directory..."
shopt -s nullglob
for ext in fa fna fasta; do
    for reference in $refdir/*/*.$ext; do
        if [[ ! -f "${reference}.bwt" ]]; then
            echo $reference >> $non_indexed_refs
        fi

        refbase=$(basename $reference)
        if [[ ! -f "$TEMPSHARE/$refbase" ]]; then
            rsync -L $reference* $TEMPSHARE
        fi
    done

    while read -r row; do
        simref=$(echo $row | col1)

        if ! grep -q $simref $non_indexed_refs; then

            if [[ ! -f "${simref}.bwt" ]]; then
                echo $simref >> $non_indexed_refs
            fi
            simref_base=$(basename $simref)
            if [[ ! -f "$TEMPSHARE/$simref_base" ]]; then
                rsync -L $simref* $TEMPSHARE
            fi

        fi
    done < $simrefs
done

#
# 3. index references
#

arr_size=$(count_lines $non_indexed_refs)

if (( arr_size != 0 )); then
    # index refs in place (copy to scratch, index there, copy index files back)
    echo "Queued: Indexing $arr_size references..."
    index_jobID=$(sbatch --array=1-$arr_size $BIN/index-reference.sh $non_indexed_refs $TEMP | col4)
else
    echo "Status: All references are already indexed, skipping indexing step"
fi

#
# 4. simulate reads
#

if [[ FLAG_SIMULATION_MODE -eq 1 ]]; then

    arr_size=$(count_lines $simrefs)

    fastq_sim_outdir="$outdir/fastq-sim"
    mkdir -p $fastq_sim_outdir

    # simulate appropriate amount of reads for each reference, work in TEMP
    echo "Queued: Simulating reads from $arr_size references"
    sim_jobID=$(sbatch --array=1-$arr_size $BIN/simulate-reads.sh $simrefs $numpairs $TEMPSHARE $TEMP | col4)

    # concat all R1s/R2s, work in TEMP, publish to fastq_sim_outdir
    concat_jobID=$(sbatch --dependency=afterok:$sim_jobID $BIN/concat-reads.sh $TEMP $fastq_sim_outdir | col4)

    # write metadata on generated reads
    python3 $BIN/write-metadata.py $simrefs $numpairs $fastq_sim_outdir

    echo "Status: Awaiting read simulation completion..."
    await_job_completion $concat_jobID
    echo "Status: Finished read simulation"
    
    read1=$TEMP/merged_raw_1.fastq
    read2=$TEMP/merged_raw_2.fastq
    # already defined if not in simulation mode
fi

#
# 5. filter/align against univec
#

# align/filter raw reads against univec, work in TEMP, copy results to fastq_sim_outdir
echo "Queued: Filtering reads against univec..."
filter_univec_jobID=$(sbatch --dependency=afterok:$concat_jobID $BIN/filter-against-univec.sh $read1 $read2 $TEMP $fastq_sim_outdir | col4)

filtered1=$fastq_sim_outdir/filtered_1.fastq
filtered2=$fastq_sim_outdir/filtered_2.fastq

#
# 5. align filtered reads to references
#

reflist=$TEMP/reflist.txt
touch $reflist

for ext in fa fna fasta; do
    for ref in $refdir/*/*.$ext; do
        echo $ref >> $reflist
    done
done

arr_size=$(count_lines $reflist)
bamdir="$outdir/bamfiles"
mkdir -p $bamdir


echo "Queued: Aligning reads to all references..."
# align filtered reads to each reference in parallel, work in TEMP for speed, publish outputs to bamdir
if [[ -s $non_indexed_refs ]]; then
    align_jobID=$(sbatch --dependency=afterok:$index_jobID:$filter_univec_jobID --array=1-$arr_size $BIN/align-reads.sh $reflist $filtered1 $filtered2 $TEMP $TEMPSHARE $bamdir | col4)
else
    align_jobID=$(sbatch --dependency=afterok:$filter_univec_jobID --array=1-$arr_size $BIN/align-reads.sh $reflist $filtered1 $filtered2 $TEMP $TEMPSHARE $bamdir | col4)
fi

#
# 6. extract top alignments
#

analysis_dir="$outdir/analysis"
mkdir -p $analysis_dir

# assign reads, work in TEMP for speed, publish outputs to analysis
echo "Queued: Assigning reads to reference of origin by alignment score..."
assign_jobID=$(sbatch --dependency=afterok:$align_jobID $BIN/assign-reads.sh $BIN/assign-reads.py $TEMP $analysis_dir | col4)

await_job_completion $assign_jobID

echo "Status: Alignment + assignment completed"
echo "done :)"
