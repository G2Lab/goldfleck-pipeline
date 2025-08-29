#!/usr/bin/env bash

#SBATCH --cpus-per-task=4
#SBATCH --mem=12G
#SBATCH --kill-on-invalid-dep=yes
#SBATCH --output=logs/slurm-assign-%A.out


script=$1
TEMP=$2
outdir=$3

module load python/3.12.3-GCCcore-13.3.0

if ! python -c "import polars" &> /dev/null; then
    pip install polars
fi

# get list of bams made in temp

touch $TEMP/bamlist.txt
for bam in $TEMP/bams/*; do
    echo $bam >> $TEMP/bamlist.txt
done

mkdir -p $TEMP/analysis

python $script -i $TEMP/bamlist.txt -o $TEMP/analysis

rsync -a $TEMP/analysis/* $outdir
