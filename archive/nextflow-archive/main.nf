#!/usr/bin/env nextflow

nextflow.enable.dsl=2


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


process WRITE_READS_METADATA {
    publishDir "$params.fastq_dir/metadata.tsv", mode: "copy"
    
    input:
    tuple path(filepath), val(num_pairs)

    output:
    path "metadata.tsv"

    script:
    """
    #!/usr/bin/env python3

    outlines = ['bamfile_path\tpercent_of_total_reads\tnumber_of_pairs_generated']

    with open("$filepath", "r") as fp:
        lines = fp.read().strip().split('\n')
        delim = ',' if ',' in lines[0] else '\t'

        for line in lines:
            fp, percent = line.split(delim)
            fp, percent, numpairs = fp, int(percent), int(percent) * ($num_pairs // 100)
            outlines.append(f"{fp}\t{percent}\t{numpairs}")

    with open("metadata.tsv", "w") as fp:
        fp.write('\n'.join(outlines) + '\n')
    """
}


process DOWNLOAD_UNIVEC {
    tag "Download UniVec_Core"

    output:
    path "UniVec_Core.fa"

    script:
    """
    univec_dir="${params.data}/univec"
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


process INDEX_REFERENCE {

    tag "$reference"

    input:
    path reference

    output:
    tuple path(reference), path("${reference}.amb"), path("${reference}.ann"), path("${reference}.bwt"), path("${reference}.pac"), path("${reference}.sa")

    script:
    """
    module load bwa/0.7.18-GCCcore-12.3.0
    bwa index $reference
    """
}



process ALIGN_TO_REFERENCE {

    tag "${reference.getBaseName()}.bam"
    publishDir "$params.bam_dir", mode: "copy"

    input:
    tuple path(fastq1), path(fastq2), path(reference), path(amb), path(ann), path(bwt), path(pac), path(sa)

    output:
    path "${reference.getBaseName()}.bam"

    script:
    """
    module load bwa/0.7.18-GCCcore-12.3.0
    module load samtools/1.21

    bwa mem -t 16 -M $reference $fastq1 $fastq2 | samtools view -bSho "${reference.getBaseName()}.bam" -
    """
}


process EXTRACT_TOP_ALIGNMENTS {
    publishDir "$params.analysis_dir", mode: 'copy'

    input:
    path bamlist

    output:
    path "*.csv"

    script:
    """
    module load python/3.12.3-GCCcore-13.3.0

    python3 ${params.bin}/assign_reads.py -i $bamlist -o "./"
    """
}


def getLocalTimestamp() {
    return java.time.LocalDateTime.now().format(
        java.time.format.DateTimeFormatter.ofPattern("yyyy-MM-dd_HH-mm-ss")
    )
}

// arguments, overwritten at runtime
params.refdir = null
params.sim_refs = null
params.num_pairs = 1_000_000
params.R1 = null
params.R2 = null
params.run_name = null

// global internal vars
params.results = params.run_name ? "$projectDir/results/" + getLocalTimestamp() + "-$params.run_name" : "$projectDir/results/" + getLocalTimestamp()
params.bin = "$projectDir/bin"
params.fastq_dir = "${params.results}/fastq-sim"
params.bam_dir = "${params.results}/bamfiles"
params.analysis_dir = "${params.results}/analysis"
params.univec_fasta = "data/univec/UniVec_Core.fa"

workflow {

    /* 1. validate arguments */

    def use_passed_fastqs = params.R1 && params.R2 && file(params.R1).exists() && file(params.R2).exists()
    def generate_fastqs   = params.sim_refs && file(params.sim_refs).exists()

    if ( ! params.refdir || ! file(params.refdir).exists() ) {
        error "Error: You must provide a reference directory with FASTA files to align to."
    }

    if ( params.univec_fasta && ! file(params.univec_fasta).exists() ) {
        error "Error: You must provide a valid filepath to --univec_fasta."
    }

    /* 2. extract & index references (channel factory used to deal with symlinked inputs) */

    Channel
        .fromPath( params.refdir )
        .flatMap { dir -> 
            def species_files = []
            def fasta_exts = ['.fa', '.fna', '.fasta']
            dir.eachFileRecurse { fp ->
                if ( fasta_exts.any { ext -> fp.name.endsWith(ext) } ) {
                    species_files.add( fp )
                }
            }
            return species_files
        }
        | INDEX_REFERENCE
        
    /* 3. Get input fastqs (whether passed or generated) */

    if ( use_passed_fastqs ) {
        Channel
            .of([file(params.R1), file(params.R2)])
            .set { fastq_pair_ch }
    }
    else if ( generate_fastqs ) {
        def delimiter = params.sim_refs.endsWith("csv") ? ',' : '\t'

        Channel
            .fromPath( params.sim_refs )
            .splitCsv( sep: delimiter )
            .map { row -> tuple( file(row[0]), row[1] * (params.num_pairs / 100) ) }
            | SIM_READS

        Channel
            .fromPath( params.sim_refs )
            .splitCsv( sep: delimiter )
            .map { row -> tuple( row[0], row[1], )}

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

    /* 4. Filter everything against & index univec */

    if ( ! file(params.univec_fasta).exists() ) {
        DOWNLOAD_UNIVEC()
    }    

    Channel
        .fromPath(params.univec_fasta)
        .concat( fastq_pair_ch )
        .collect()
        | FILTER_AGAINST_UNIVEC

    /* 4. align reads */

    FILTER_AGAINST_UNIVEC.out
        .combine( INDEX_REFERENCE.out )
        | ALIGN_TO_REFERENCE

    /* 5. extract alignment data from bamfiles */

    ALIGN_TO_REFERENCE.out
        .collectFile( name: 'bamlist.txt', newLine: true )
        | EXTRACT_TOP_ALIGNMENTS
}   