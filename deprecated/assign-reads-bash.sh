#!/bin/bash

#SBATCH --job-name=assign-reads
#SBATCH --time=144:00:00
#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --mem=30G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=email

module load samtools

#get the good reads
dir=/gpfs/commons/home/apandit/gghome/lions/align-out-mut0/bamfiles

a=($(ls ${dir}/*.bam))
b=($(ls ${dir}/*.bam | wc -l))


# samtools view -F 260 $dir/lion.bam | awk '{if ($5>=20) {print $1"\t"$5"\t"$12}}' | awk -F":|\t" '{if ($8>=20 && $8<30 && $11<6) {print $1":"$2":"$3":"$4":"$5":"$6":"$7}}' | awk '!seen[$0]++' >  $dir/lion.bam.1.txt
# samtools view -F 260 $dir/lion.bam | awk '{if ($5>=20) {print $1"\t"$5"\t"$12}}' | awk -F":|\t" '{if ($8>30 && $11<4) {print $1":"$2":"$3":"$4":"$5":"$6":"$7}}' | awk '!seen[$0]++' >>  $dir/lion.bam.2.txt


for ((i=0; i<$b; i++)); do
	samtools view -F 260 ${a[$i]} | awk '{if ($5>=20) {print $1"\t"$5"\t"$12}}' | awk -F":|\t" '{if ($8>=20 && $8<30 && $11<6) {print $1":"$2":"$3":"$4":"$5":"$6":"$7}}' | awk '!seen[$0]++' >  ${a[$i]}.txt
	samtools view -F 260 ${a[$i]} | awk '{if ($5>=20) {print $1"\t"$5"\t"$12}}' | awk -F":|\t" '{if ($8>=30 && $11<4) {print $1":"$2":"$3":"$4":"$5":"$6":"$7}}' | awk '!seen[$0]++' >>  ${a[$i]}.txt
done

for ((i=0; i<$b; i++)); do
	k=$(( $i+1 ))
	for ((j=$k; j<$b; j++)); do
		awk 'NR==FNR{seen[$0]=1; next} seen[$0]' ${a[$i]}.txt ${a[$j]}.txt > ${dir}/tmps
		if [[ $(wc -l <${dir}/tmps) -ge 1 ]]; then
			awk 'NR==FNR{seen[$0]=1; next} !seen[$0]' ${a[$i]}.txt ${a[$j]}.txt > tmp1
			mv tmp1 ${a[$j]}.txt
			awk 'NR==FNR{seen[$0]=1; next} !seen[$0]' ${a[$j]}.txt ${a[$i]}.txt > tmp1
			mv tmp1 ${a[$i]}.txt
		fi
	done
done