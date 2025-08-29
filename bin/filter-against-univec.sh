#!/usr/bin/env bash

#SBATCH --cpus-per-task=2
#SBATCH --mem=4G
#SBATCH --output=logs/slurm-filterunivec-%A.out


merged_1=$1
merged_2=$2
TEMP=$3
fastq_dir=$4

if [[ -z $TEMP ]]; then
    exit 1
fi

module load picard/3.0.0-Java-17
module load samtools/1.21
module load bwa/0.7.18-GCCcore-12.3.0

picard=$EBROOTPICARD/picard.jar
univec=$TEMP/UniVec_Core.fa

bwa index $univec
bwa aln -n 1 -t 4 $univec $merged_1 > $TEMP/samplePair1.sai
bwa aln -n 1 -t 4 $univec $merged_2 > $TEMP/samplePair2.sai
bwa sampe $univec $TEMP/samplePair1.sai $TEMP/samplePair2.sai $merged_1 $merged_2 | samtools view -bS > $TEMP/univec_aligned.bam 

samtools view -h -f 4 $TEMP/univec_aligned.bam | col1 > $TEMP/univec_unmapped.txt
java -Xmx200g -jar $picard FilterSamReads -VALIDATION_STRINGENCY LENIENT -I $TEMP/univec_aligned.bam -O $TEMP/output.bam -READ_LIST_FILE $TEMP/univec_unmapped.txt -FILTER includeReadList

samtools sort -n $TEMP/output.bam -o - | samtools fastq -1 $TEMP/filtered_1.fastq -2 $TEMP/filtered_2.fastq

echo "DEBUG: rsync $TEMP/filtered_*.fastq $fastq_dir"
rsync $TEMP/filtered_*.fastq $fastq_dir
