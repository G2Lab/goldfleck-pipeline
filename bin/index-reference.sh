#!/usr/bin/env bash

#SBATCH --cpus-per-task=2
#SBATCH --mem=10G
#SBATCH --output=logs/slurm-index-%A_%a.out


non_indexed_refs=$1
TEMP=$2

refpath=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$non_indexed_refs")
refname=$(basename $refpath)
refdir=$(dirname $refpath)

cp $refpath $TEMP/$refname

module load bwa/0.7.18-GCCcore-12.3.0
bwa index $TEMP/$refname

cp $TEMP/$refname.* $refdir
