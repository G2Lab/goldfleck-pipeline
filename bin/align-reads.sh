#!/usr/bin/env bash

#SBATCH --cpus-per-task=16
#SBATCH --mem=48G
#SBATCH --kill-on-invalid-dep=yes
#SBATCH --output=logs/slurm-align-%A_%a.out

reflist=$1
fastq1=$2
fastq2=$3
TEMP=$4
TEMPSHARE=$5
outdir=$6

if [[ "$#" -ne 6 ]]; then
    echo "incorrect argument count"
    exit 1
fi


refpath=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$reflist")
refname=$(basename $refpath)
refpath=$TEMPSHARE/$refname

refstem="${refname%.*}"

echo "DEBUG"
echo "ref=$ref"
echo "fastq1=$fastq1"
echo "fastq2=$fastq2"
echo "refstem=$refstem"

module load bwa/0.7.18-GCCcore-12.3.0
module load samtools/1.21

mkdir -p $TEMP/bams

bwa mem -t 16 -M $refpath $fastq1 $fastq2 | samtools view -bSho $TEMP/bams/$refstem.bam -

rsync $TEMP/bams/$refstem.bam $outdir
