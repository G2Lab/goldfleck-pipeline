#!/usr/bin/bash

# run filter-and-align.sh on simulated lion fastqs

script="/gpfs/commons/home/apandit/gghome/lions/scripts/filter-and-align.sh"
refdirs="/gpfs/commons/home/apandit/gghome/lions/genomes"
univec="/gpfs/commons/home/apandit/gghome/lions/univec/UniVec_Core"

# below: already ran w/ mutrate=0, include if not
for mutrate in 0 0.001 0.0025 0.005 0.01 0.025 0.05; do
    fastq1="/gpfs/commons/home/apandit/gghome/lions/wgsim-out/fastqs/lionref/lionref-mut${mutrate}-1.fq"
    fastq2="/gpfs/commons/home/apandit/gghome/lions/wgsim-out/fastqs/lionref/lionref-mut${mutrate}-2.fq"

    echo "[run-alignment] Queueing filter-and-align.sh on the cluster with paired fastqs lionref-mut${mutrate}"
    sbatch $script -1 $fastq1 -2 $fastq2 -d $refdirs -u $univec -o align-out-mut$mutrate
done
