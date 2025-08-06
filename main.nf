#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// workaround to params.refdir/*/*.fa not capturing symlink *.fa
process COLLECT_REFERENCES {
	output:
	path "references.txt"

	script:
	"""
	shopt -s nullglob
	echo -n "" > references.txt
	for ext in fa fna fasta; do
		for file in ${file(params.refdir)}/*/*.\$ext; do
			echo \$file >> references.txt
		done
	done
	"""
}


process SIM_READS {
    input:
    tuple path(fasta), val(num_pairs)

    output:
    tuple path("${fasta.baseName}_1.fastq"), path("${fasta.baseName}_2.fastq")

    script:
    """
    module load samtools/1.21

    file_name="$fasta"
    base_name="\${file_name%.*}"
    wgsim $fasta \${base_name}_1.fastq \${base_name}_2.fastq -e 0.02 -d 100 -s 30 -N $num_pairs -1 150 -2 150 -r 0.001
    """
}


process CONCAT_READS {
    publishDir "$params.fastq_dir", mode: "copy"

    input:
    tuple val(Rtype), path(files)

    output:
    path "merged_raw_${Rtype}.fastq"

    script:
    """
    cat ${files.join(' ')} > merged_raw_${Rtype}.fastq
    """
}


process DOWNLOAD_UNIVEC {
    tag "Download UniVec_Core"

    output:
    path "UniVec_Core.fa"

    script:
    """
    univec_dir="$projectDir/data/univec"
    mkdir -p \$univec_dir
    curl -L -o \$univec_dir/UniVec_Core.fa ftp://ftp.ncbi.nlm.nih.gov/pub/UniVec/UniVec_Core
    ln -s \$univec_dir/UniVec_Core.fa UniVec_Core.fa
    """
}


process FILTER_AGAINST_UNIVEC {
    publishDir "$params.fastq_dir", mode: "copy"

    input:
    tuple path(univec), path(merged_1), path(merged_2)

    output:
    tuple path("filtered_1.fastq"), path("filtered_2.fastq")

    script:
    """
    module load picard/3.0.0-Java-17
    module load samtools/1.21
    module load bwa/0.7.18-GCCcore-12.3.0
    
    picard=\$EBROOTPICARD/picard.jar

	if [[ ! -f ${univec}.bwt ]]; then
		bwa index $univec
	fi

    tmp="/nfs/scratch/\$SLURM_JOB_ID"
    mkdir -p \$tmp
    bwa aln -n 1 -t 4 $univec $merged_1 > \$tmp/samplePair1.sai
    bwa aln -n 1 -t 4 $univec $merged_2 > \$tmp/samplePair2.sai
    bwa sampe $univec \$tmp/samplePair1.sai \$tmp/samplePair2.sai $merged_1 $merged_2 | samtools view -bS > \$tmp/univec_aligned.bam 

    samtools view -h -f 4 \$tmp/univec_aligned.bam | awk '{print \$1}' > \$tmp/univec_unmapped.txt
    java -Xmx200g -jar \$picard FilterSamReads -VALIDATION_STRINGENCY LENIENT -I \$tmp/univec_aligned.bam -O \$tmp/output.bam -READ_LIST_FILE \$tmp/univec_unmapped.txt -FILTER includeReadList

    samtools sort -n \$tmp/output.bam -o - | samtools fastq -1 filtered_1.fastq -2 filtered_2.fastq
    """
}


process ALIGN_TO_REFERENCE {
    publishDir "$params.bam_dir", mode: "copy"

    input:
    tuple path(fastq1), path(fastq2), path(reference)

	output:
	path "${reference.getBaseName()}.bam"

	script:
	"""
	module load bwa/0.7.18-GCCcore-12.3.0
	module load samtools/1.21

	# index in the actual reference dir not the nextflow work dir
	reference=\$(readlink $reference)

    if [[ ! -f "\$reference.bwt" ]]; then
        bwa index "\$reference"
    fi

	# get a reference to all the index files in the work dir so alignment can happen there
	for ext in amb ann bwt pac sa; do
		ln -s "\$reference.\$ext" ${reference}.\$ext
	done

	bwa mem -t 16 -M $reference $fastq1 $fastq2 | samtools view -bSho "${reference.getBaseName()}.bam -"
	"""
}


// arguments, nulls are required to be passed at runtime
params.refdir = null
params.sim_refs = null
params.num_pairs = 1_000_000
params.R1 = null
params.R2 = null

// global internal vars
params.results = "$projectDir/results"
params.fastq_dir = "${params.results}/fastq-sim"
params.bam_dir = "${params.results}/bamfiles"
params.univec_fasta = "data/univec/UniVec_Core.fa"

workflow {
    // check conditionals
    def use_passed_fastqs = params.R1 && params.R2 && file(params.R1).exists() && file(params.R2).exists()
    def generate_fastqs   = params.sim_refs && file(params.sim_refs).exists()

	/* 1. get references into a channel because nextflow globbing makes me sad */

    if ( ! params.refdir || ! file(params.refdir).exists() ) {
        error "Error: You must provide a reference directory with FASTA files to align to."
    }
	COLLECT_REFERENCES()
	COLLECT_REFERENCES.out
		.splitText()
		.filter { it }
		.map { file(it.trim()) }  // remove trailing newlines
		.set { refdir_ch }

    /* 2. Get input fastqs (whether passed or generated) */

    if ( use_passed_fastqs ) {
        Channel
			.of([file(params.R1), file(params.R2)])
			.set { fastq_pair_ch }
    }
    else if ( generate_fastqs ) {
        Channel
            .fromPath( params.sim_refs )
            .splitCsv( sep: "\t" )
            .map { row -> tuple( file(row[0]), row[1] ) }
            | SIM_READS

        SIM_READS.out
            .flatten()
			.map { tuple( it.name.endsWith('1.fastq') ? '1' : '2', file(it) ) }
            .groupTuple()
            | CONCAT_READS

		CONCAT_READS.out
			.set { fastq_pair_ch }
    }
    else {
        error "Error: You must provide either --R1 and --R2 or --sim_refs (and ensure the files exist)."
    }

    /* 3. Filter everything against & index univec */

    if ( ! file(params.univec_fasta).exists() ) {
        DOWNLOAD_UNIVEC()
    }	

    Channel
		.fromPath(params.univec_fasta)
		.concat( fastq_pair_ch )
		.collect()
		| FILTER_AGAINST_UNIVEC

	/* 4. index & align against references */

	FILTER_AGAINST_UNIVEC.out
		.combine( refdir_ch )
		| ALIGN_TO_REFERENCE

	ALIGN_TO_REFERENCE.out
		.view()

	/* 5. */
}   