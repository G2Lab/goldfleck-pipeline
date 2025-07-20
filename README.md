# Multi-genome Alignment Pipeline

A pipeline designed to
- Simulate a set of paired-end reads with varying mutation rates from a set of reference genomes using `wgsim`
- Align simulated reads to each reference genome and generate `bam` files for each pair of reads, filtering out known contaminants (e.g. PCR primers) using NCBI's UniVec database
- Analyze bamfile statistics and demonstrate that bamfile from a given species `A` can be uniquely mapped back to species `A`, and not the other species

## Files

`scripts/filter-and-align.sh`

Aligns provided FASTQ files to every genome (FASTA file) provided, first discarding known contaminants in [NCBI's UniVec database](https://www.ncbi.nlm.nih.gov/tools/vecscreen/univec/). For each genome alignment, a bam file is created in the current working directory (e.g. `read1.fq`, `read2.fq`, `lion.fa` -> `lion.bam`). UniVec filtering is done in `filter-and-align.sh`, which then 

| Argument | Description |
|-------------------|-------------------------|
| `-1 <file.fastq>` | paired-end fastq file 1 |
| `-2 <file.fastq>` | paired-end fastq file 2 |
| `-d <directory>` | reference sequences (structured as directory/genome/genome.fa) |
| `-u <univec path>` | (optional) path to a univec reference file, found/downloaded if not provided |
| `-o <outdir>` | (optional) path to output directory, `align-output` if unspecified |



