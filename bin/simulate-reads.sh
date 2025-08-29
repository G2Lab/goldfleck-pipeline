#!/usr/bin/env bash

#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --output=logs/slurm-sim-%A_%a.out


file_name=$1
total_pairs=$2
TEMPSHARE=$3
TEMP=$4

# Select the line for this task
line=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$file_name")

# Extract arguments
fasta_name=$(echo "$line" | awk '{print $1}')
percent=$(echo "$line" | awk '{print $2}')

module load samtools/1.21

fasta_name=$(basename $fasta_name)
base_name="${fasta_name%.*}"
num_pairs=$(( $total_pairs * $percent / 100 ))

echo "DEBUG: wgsim $TEMPSHARE/$fasta_name $TEMP/${base_name}_1.fastq $TEMP/${base_name}_2.fastq"
echo "DEBUG: total_pairs=$total_pairs"
echo "DEBUG: num_pairs=$num_pairs"  

wgsim  $TEMPSHARE/$fasta_name $TEMP/${base_name}_1.fastq $TEMP/${base_name}_2.fastq -e 0.02 -d 100 -s 30 -N $num_pairs -1 150 -2 150 -r 0.001
