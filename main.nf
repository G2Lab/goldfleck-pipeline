#!/usr/bin/env nextflow
nextflow.enable.dsl=2

params.results = "$projectDir/results"
params.fastq_dir = "${params.results}/fastq-sim"

params.sim_refs = null
params.num_pairs = 1_000_000
params.R1 = null
params.R2 = null


process SIM_FASTQS_AND_MERGE {
    publishDir "${params.fastq_dir}", mode: 'copy', pattern: 'merged_*.fastq'

    def script_path = "/gpfs/commons/home/apandit/gghome/lions/bin/gen-fastq.py"

    input:
    path sim_refs
    val num_pairs

    output:
    path "merged_R1.fastq"
    path "merged_R2.fastq"

    script:
    """
    set -euo pipefail

    tmpdir="/nfs/scratch/\$USER/fastq-work-\$(uuidgen)"
    mkdir -p \$tmpdir

    module load python
    python $script_path -r $sim_refs -n $num_pairs -o \$tmpdir -b "nextflow-test" \\
        | grep "Submitted batch job " | awk '{print \$NF}' > \$tmpdir/job_ids.tmp

    if [[ ! -s \$tmpdir/job_ids.tmp ]]; then
        echo "Error: No job IDs captured."
        exit 1
    fi

    while read jobid; do
        echo "Waiting for job \$jobid..."
        while squeue -j \$jobid | grep -q \$jobid; do sleep 5; done
        echo "Job \$jobid done."
    done < \$tmpdir/job_ids.tmp

    cat \$tmpdir/*R1.fastq > merged_R1.fastq
    cat \$tmpdir/*R2.fastq > merged_R2.fastq

    rm -rf \$tmpdir
    """
    // flow
    // 1. set euo pipefall so no silent fails
    // 2. run the script
    // 3. validate that we have slurm processes & wait for them to finish
    // 4. merge the fastqs from each species into combined
    // 5. cleanup
}


// process ALIGN_READS {
//     // tbd, this will produce bamfiles
// }


// process ASSIGN_READS {
//     // tbd, this will assign reads and produce all those lovely csvs
// }


workflow {
    def use_real_fastqs = file(params.R1).exists() && file(params.R2).exists()
    def use_sim_refs = !use_real_fastqs && file(params.sim_refs).exists()

    // use passed fastqs or simulate if not present (with input error handling)
    if (use_real_fastqs) {
        r1_ch = Channel.fromPath(params.R1)
        r2_ch = Channel.fromPath(params.R2)
    }
    else if (use_sim_refs) {
        sim_refs_ch = Channel.fromPath(params.sim_refs)
        num_pairs_ch = Channel.of(params.num_pairs)

        sim_out = SIM_FASTQS_AND_MERGE(sim_refs_ch, num_pairs_ch)

        r1_ch = sim_out.out[0]
        r2_ch = sim_out.out[1]
    }
    else {
        error "Error: please either pass --R1 and --R2 or --sim_refs. If you did, ensure the files exist (you can run `file {filepath}` to check)"
    }


}
