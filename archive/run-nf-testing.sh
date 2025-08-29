#!/bin/bash

PROJ="$HOME/gghome/lions"

module load nextflow

if [[ -f "all_hashes.txt" ]]; then rm "all_hashes.txt"; fi
touch "all_hashes.txt"

for blend in blend1 blend2 blend3; do 
    nextflow run $PROJ/main.nf --refdir $PROJ/testing/genomes/ --sim_refs $PROJ/testing/$blend.tsv --univec_fasta $PROJ/data/univec/UniVec_Core.fa --run_name "$blend" -dump-hashes >> all_hashes.txt
done
