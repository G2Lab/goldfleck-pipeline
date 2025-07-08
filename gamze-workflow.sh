#!/bin/bash

#SBATCH -p cpu
#SBATCH -J workflow
#SBATCH -n 1 --mem 300000 -t 48:00:00 
#SBATCH --mail-type=ALL
#SBATCH --mail-user=email


# This is the script to map the raw sequences 
# to the database of reference genomes and mark
# PCR duplicates.
# This script creates N number of jobs to be
# submitted to the server to run every species
# in parallel.
# by Gamze Gursoy
# edited by Akash Pandit
# last edit: July 2nd, 2025

# load modules
module load picard/3.0.0-Java-17
module load samtools/1.21



# Step 1: Find samples to process
sampleDir=$1

echo "testing"

declare -a fqPaths

shopt -s nullglob

i=0
for fp in "${sampleDir}"/*; do
    for ext in fastq.gz fastq fq.gz fq; do
        if [[ $fp == *.$ext ]]; then
            fqPaths[$i]="$fp"
            ((i++))
            break
        fi
    done
done


if [[ ${i} == 0 ]]; then
	echo "No fastq file has been found in this folder, please input the correct folder"
	exit 1
fi

echo "pair1 is found: "${fqPaths[0]}
echo "pair2 is found: "${fqPaths[1]}

# Step 2: Filter out known contaminant reads by mapping to Univec

# specify tools/directories
picard=/ysm-gpfs/project/miranker/tools/picard_tools/picard.jar
bwa=/ysm-gpfs/project/miranker/tools/bwa-0.7.15/bwa
univecdir=/ysm-gpfs/project/miranker/workflow/univec/
univecindex=UniVecIndex
tmp=/ysm-gpfs/scratch60/gg487/Sample_MRG-191217-1_189_361_05March/tmp/

# align to univec & ID non-aligned reads to keep
mkdir ${tmp}
${bwa} aln -n 1 ${univecdir}/${univecindex} ${fqPaths[0]} > ${tmp}/saiFile1.sai
${bwa} aln -n 1 ${univecdir}/${univecindex} ${fqPaths[1]} > ${tmp}/saiFile2.sai
${bwa} sampe ${univecdir}/${univecindex} ${tmp}/saiFile1.sai ${tmp}/saiFile2.sai ${fqPaths[0]} ${fqPaths[1]} | samtools view -bS > ${tmp}/mapped.bam
samtools view -h -f 4 ${tmp}/mapped.bam | awk '{print $1}' > ${tmp}/unmapped.txt
java -Xmx200g -jar ${picard} FilterSamReads VALIDATION_STRINGENCY=LENIENT I=${tmp}/mapped.bam O=${tmp}/output.bam READ_LIST_FILE=${tmp}/unmapped.txt FILTER=includeReadList  
samtools sort -n ${tmp}/output.bam -o ${tmp}/output.s.bam
samtools fastq -1 ${tmp}/paired1.fq -2 ${tmp}/paired2.fq -n ${tmp}/output.s.bam
rm -rf tmp

# point to fastqs that have passed above filtering step
fqPaths[0]=${tmp}/paired1.fq
fqPaths[1]=${tmp}/paired2.fq

echo "pair1 is found: "${fqPaths[0]}
echo "pair2 is found: "${fqPaths[1]}

# Step 4: Create the jobs that map and mark duplicates

# define output file
output=${sampleDir}/BAMs
if [ ! -d $output ]; then
	mkdir $output
fi


#Define the directory to the database as the second user argument
dir2=$2

#iterate over N reference_genomes in the database and submit jobs for each mapping
declare -a reference_genomes
N=0
for f in ${dir2}/*; do
    if [ -d ${f} ]; then
		reference_genomes[$N]=${f}
		N=$(( $N+1 ))    
    fi
done

cd ${output}

for ((i=0; i<$N; i++))
do
	echo "#!/bin/sh
#SBATCH -p cpu
#SBATCH -J bam_$i -t 240:00:00
#SBATCH -c 16 --mem-per-cpu=10000

	module load SAMtools
        samtools=\$EBROOTSAMTOOLS/bin/samtools

	#find the reference file
	for i in ${reference_genomes[$i]}/*; do
        if [[ \${i} == *.fna.gz ]];then
			gunzip \${i}
		elif [[ \${i} == *.fasta.gz ]];then
            gunzip \${i}
		fi
	done

	for i in ${reference_genomes[$i]}/*; do  # gets last .fna or .fasta file path
		if [[ \${i} == *.fna ]];then
			ref=\${i}
		elif [[ \${i} == *.fasta ]];then
            ref=\${i}
		fi
	done

	#check if the reference files are indexed
 	indicator=0
	for i in ${reference_genomes[$i]}/*; do
    		if [[ \${i} == *.bwt ]];then
			indicator=0
			break
		else
			indicator=1
		fi
	done

	if [  \${indicator} == 1 ];then
		${bwa} index \$ref
	fi

	#get the name of the species
	IFS='/' read -ra arr <<< "${reference_genomes[$i]}"
	len=\${#arr[@]}
	species=\${arr[\$(( \${len}-1 ))]}

	#map the fastqs to reference genome
	${bwa} mem -t 16 -M \$ref ${fqPaths[0]} ${fqPaths[1]} | \${samtools} view -bSho \${species}.bam -" > $i.sh
	sbatch $i.sh 
done 

assumed structure

reference_genomes/
├── species1/
│	├── species1.fasta
│	├── species1.fasta.bwt
│	...
├── species2/
...	├── species2.fasta
	...
