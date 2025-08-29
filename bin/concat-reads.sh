#!/usr/bin/env bash

#SBATCH --cpus-per-task=1
#SBATCH --mem=1G
#SBATCH --kill-on-invalid-dep=yes
#SBATCH --output=logs/slurm-concat-%A.out


TEMP=$1
fastq_dir=$2

if [[ "$#" -lt 2 ]]; then
    echo "missing arguments"
    exit 1
fi

shopt -s nullglob
cat $TEMP/*1.fastq > $TEMP/merged_raw_1.fastq
cat $TEMP/*2.fastq > $TEMP/merged_raw_2.fastq

if [[ ! -s $TEMP/merged_raw_1.fastq || ! -s $TEMP/merged_raw_2.fastq ]]; then
    echo "Error: a merged file has no reads, something's fishy here"
    exit 1
fi

echo "DEBUG: rsync $TEMP/merged_raw_*.fastq $fastq_dir"
rsync "$TEMP"/merged_raw_*.fastq $fastq_dir
