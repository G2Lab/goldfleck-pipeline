#!/usr/bin/bash

#SBATCH --job-name=filter-univec
#SBATCH --time=48:00:00
#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --mem=20G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=email

# Filter paired end fastq files against univec contaminant
# database (mark pcr duplicates) and align remaining reads
# to all reference genomes found in the passed directory
#
# Creates N jobs for N references that run in parallel
# By Akash Pandit
# Last edited: July 3rd, 2025
#
# Workflow:
# 1. parse command line arguments
# 2. load required modules
# 3. map to univec
# 4. generate & run slurm scripts mapping to each ref seq

####################################################
# 1. Parse command line arguments
####################################################

# Usage message
usage() {
    echo "Usage: $0 -1 <read1.fastq.gz> -2 <read2.fastq.gz> -d <refdirs> [-u <univec>]"
    echo ""
    echo "  -1  Path to R1 FASTQ file"
    echo "  -2  Path to R2 FASTQ file"
    echo "  -d  Directory containing reference directories (each containing the ref FASTA file)"
    echo "  -u  Path to UniVec (or UniVec Core) fasta file, if present"
    echo "  -h  Show this help message"
    exit 1
}
[[ $# -eq 0 ]] && usage  # display usage on no flag run

# Initialize variables
read1=""
read2=""
refdirs=""

# Parse flags
while getopts ":1:2:d:u:h" opt; do
    case $opt in
        1) read1="$OPTARG" ;;
        2) read2="$OPTARG" ;;
        d) refdirs="$OPTARG" ;;
        u) univec="$OPTARG" ;;
        h) usage ;;
        \?) echo "Error: Invalid option -$OPTARG" >&2; usage ;;
        :) echo "Error: Option -$OPTARG requires an argument." >&2; usage ;;
    esac
done

# Validate required arguments
if [[ -z "$read1" || -z "$read2" || -z "$refdirs" ]]; then
    echo "Error: Missing required arguments." >&2
    usage
elif [[ ! -f "$read1" ]]; then
    echo -e "\nError: file $read1 does not exist.\n" >&2
    exit 1
elif [[ ! -f "$read2" ]]; then
    echo -e "\nError: file $read2 does not exist.\n" >&2
    exit 1
elif [[ ! -d "$refdirs" ]]; then
    echo -e "\nError: path $refdirs is not a directory.\n" >&2
    exit 1
elif [[ -n "$univec" && ! -f "$univec" ]]; then
    echo -e "\nError: invalid univec path passed. File $univec does not exist.\n" >&2
    exit 1
fi
refdirs=$(realpath $refdirs)  # 

####################################################
# 2. load required modules & other dependencies
####################################################

module load picard/3.0.0-Java-17
module load samtools/1.21
module load bwa/0.7.18-GCCcore-12.3.0

picard=$EBROOTPICARD/picard.jar

if [[ ! -f "$univec" ]]; then
    univec=$(realpath UniVec_Core)
    echo "Status: UniVec_Core not found in current working directory, downloading from ncbi..."
    curl https://ftp.ncbi.nlm.nih.gov/pub/UniVec/UniVec_Core > $univec
    if [ $? -ne 0 ]; then
        echo "Status: could not reach ncbi, please check your internet connection"
        exit 1
    fi
fi

####################################################
# 3. align to univec
####################################################

if [[ ! -f "$univec.bwt" ]]; then
    bwa index $univec   # -> cwd/*
fi

tmp="tmp-$SLURM_JOB_ID"
mkdir $tmp
bwa aln -n 1 -t 4 $univec $read1 > $tmp/samplePair1.sai  # -> cwd/tmp/samplePair1.sai
bwa aln -n 1 -t 4 $univec $read2 > $tmp/samplePair2.sai  # -> cwd/samplePair2.sai
bwa sampe $univec $tmp/samplePair1.sai $tmp/samplePair2.sai $read1 $read2 | samtools view -bS > $tmp/univec_aligned.bam    # cwd/univec_aligned.bam

samtools view -h -f 4 $tmp/univec_aligned.bam | awk '{print $1}' > $tmp/univec_unmapped.txt
java -Xmx200g -jar $picard FilterSamReads -VALIDATION_STRINGENCY LENIENT -I $tmp/univec_aligned.bam -O $tmp/output.bam -READ_LIST_FILE $tmp/univec_unmapped.txt -FILTER includeReadList

samtools sort -n $tmp/output.bam -o $tmp/output.s.bam
mkdir fastqs

read1=$(basename $read1)
read2=$(basename $read2)
paired1="${read1%.*}-paired1.fq"
paired2="${read2%.*}-paired2.fq"
read1=$(realpath fastqs/$paired1)
read2=$(realpath fastqs/$paired2)

samtools fastq -n $tmp/output.s.bam -1 $read1 -2 $read2
rm -rf tmp



####################################################
# 4. create & run slurm scripts
####################################################

if [[ ! -d "./slurm" ]]; then
    mkdir ./slurm
fi


for refdir in $refdirs/*; do
    refdir_name=$(basename $refdir)
    scriptname=slurm/slurm-align-$refdir_name.sh

    echo "#!/bin/bash
#SBATCH -p gglab_cpu 
#SBATCH -J bam_$refdir_name -t 240:00:00
#SBATCH -c 16 --mem-per-cpu=10000
#SBATCH -o slurm/slurm-%j.out

module load samtools/1.21

#find the reference file

for f in $refdir/*; do
    if [[ \"\$f\" == *.fna.gz || \"\$f\" == *.fasta.gz || \"\$f\" == *.fa.gz ]]; then
        gunzip \$f
    fi
done

for f in $refdir/*; do  # gets last .fna or .fasta file path
    if [[ \"\$f\" == *.fna || \"\$f\" == *.fasta || \"\$f\" == *.fa ]]; then
        ref=\${f}
    fi
done

#check if the reference files are indexed

no_index_found=1
for f in $refdir/*; do
    if [[ \"\$f\" == *.bwt ]];then
        no_index_found=0
    fi
done

if [  \${no_index_found} -eq 1 ];then
    echo \$ref
    bwa index \$ref
fi

#map the fastqs to reference genome
bwa mem -t 16 -M \$ref ${read1} ${read2} | samtools view -bSho ${refdir_name}.bam -" > $scriptname
    sbatch $scriptname
done 