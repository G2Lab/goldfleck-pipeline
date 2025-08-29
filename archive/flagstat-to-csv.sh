#!/bin/bash

# Parse the output of samtools flagstat (only QC-passed reads)
# of multiple bam files into csv for ease of plotting
#
# By: Akash Pandit
# Last Modified: 07/15/2025

module load samtools

for mutrate in 0 0.001 0.0025 0.005 0.01 0.025 0.05; do
    bampath="/gpfs/commons/home/apandit/gghome/lions/align-out-mut$mutrate/bamfiles"
    csv_passed="pass-qc/flagstat_qcpass_$mutrate.csv"
    csv_failed="fail-qc/flagstat_qcfail_$mutrate.csv"

    # setup column names
    for csvpath in $csv_passed $csv_failed; do
        echo -n "file_path,"        >  $csvpath
        echo -n "total_reads,"      >> $csvpath
        echo -n "primary,"          >> $csvpath
        echo -n "secondary,"        >> $csvpath
        echo -n "suppelmentary,"    >> $csvpath
        echo -n "dupes,"            >> $csvpath
        echo -n "primary_dupes,"    >> $csvpath
        echo -n "mapped,"           >> $csvpath
        echo -n "primary_mapped,"   >> $csvpath
        echo -n "paired_in_seq,"    >> $csvpath
        echo -n "read1count,"       >> $csvpath
        echo -n "read2count,"       >> $csvpath
        echo -n "properly_paired,"  >> $csvpath
        echo -n "self_mate_mapped," >> $csvpath
        echo -n "singletons,"       >> $csvpath
        echo -n "mate_mapped_diffchr_mapqLT5," >> $csvpath
        echo    "mate_mapped_diffchr_mapqGE5"  >> $csvpath
    done

    # populate with data
    for file in $bampath/*.bam; do
        echo -n $file >> $csv_passed
        echo -n $file >> $csv_failed

        samtools flagstat $file > flagstat.tmp
        cat flagstat.tmp | awk '{printf ",%s", $1} END {print ""}' >> $csv_passed
        cat flagstat.tmp | awk '{printf ",%s", $3} END {print ""}' >> $csv_failed
        rm flagstat.tmp
        echo "[Status] finished $file"
    done
done