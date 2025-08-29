# Multi-genome Alignment Pipeline (Goldfleck Pipeline)

A pipeline designed to:

- (Optionally) Simulate a set of paired-end reads with varying mutation rates from a set of reference genomes using `wgsim`
- Align simulated reads to each reference genome and generate `bam` files for each pair of reads, filtering out known 
  contaminants (e.g. PCR primers) using NCBI's UniVec database
- Analyze bamfile statistics and demonstrate that bamfile from a given species `A` can be uniquely mapped back to species `A`, 
  and not the other species

## Usage

This pipeline can either be ran with a set of post-QC paired-end reads or with a TSV containing a 'recipe' for simulating FASTQ
generation with `wgsim`. Note that arguments in square brackets are optional. Running this 

Running with existing FASTQs: `bash main.nf  -d <reference directory> -1 <first paired FASTQ> -2 <second paired FASTQ>`

Running with FASTQ simulation: `bash main.nf -d <reference directory> -s <path to CSV/TSV recipe> [-n <read pair count>]`

## Arguments

| Argument | Description | Optional | Default Value |
|-----|-----|-----|-----|
| `-d <path>` | Path to reference directory, structured as `refdir/genome/genome.{fa,fna,fasta}` | No | N/A |
| `-1 <path>` | Path to first paired FASTQ file | Not without simulation | N/A |
| `-2 <path>` | Path to second paired FASTQ file | Not without simulation | N/A |
| `-s <path>` | Path to CSV or TSV specifying simulated FASTQ composition | Not without -1/-2 | N/A |
| `-n <int>` | Number of read pairs to simulate, only affects output for simulation | Yes | 1,000,000 |
| `-r <run label>` | Label to be appended to output directory, `results/timestamp` -> `results/timestamp-label` | Yes | N/A |
| `-h` | Print usage instructions & argument descriptions, then exit | Yes | N/A |

#### Note on `-s`:
The CSV/TSV file must contain 2 columns and no header, with column 1 containing paths to a reference FASTA, and column 2 containing 
percentage 'weights' for how much each reference will contribute to the overall FASTQ pair. As these weight values are percents,
they must be integer values which add up to 100, all within the 0 to 100 range.  **Paths must be absolute, relative paths may cause
undefined behavior.**

## Pipeline Output

<!-- box drawing characters: ─ └ ├ │ components └── ├── │   -->
```
results/
├── timestamped-run-output/
│   │
│   ├── analysis/
│   │   ├── species1-unique-raw.csv
│   │   ├── species1-unique-top-alignments.csv
│   │   ├── ...
│   │   └── top-alignments.csv
│   │   
│   ├── bamfiles/
│   │   ├── species1.bam
│   │   └── ...
│   │
│   └── fastq-sim/
│       ├── blend_metadata.tsv
│       ├── filtered_1.fastq
│       ├── filtered_2.fastq
│       ├── merged_raw_1.fastq
│       └── merged_raw_2.fastq
│    
...
```

### `timestamped-run-output/`
Each run of the pipeline will generate its results in its own timestamped directory, formatted as `YYYY-MM-DD_hh-mm-ss`.

### `fastq-sim/`
This subdirectory contains FASTQ files after univec filtering and if applicable, simulation. Simulated fastqs are 
`merged_raw_{1,2}.fastq`, with post-univec filtering being `filtered_{1,2}.fastq`. 

If simulated fastqs are generated, a `blend_metadata.tsv` file will be generated as well, detailing what proportion and 
the raw number of read pairs each reference fasta contributed to the pair.

### `bamfiles/`
This subdirectory contains the bam files generated after aligning `filtered_{1,2}.fastq` in the above directory to each
reference genome provided. The bamfiles inherit the base file name of their reference files.

### `analysis/`
This subdirectory contains the analysis results of each bamfile. A global file top-alignments.csv is generated, where reads
that met MAPQ and NM quality criteria for *all references* are included, with MAPQ, NM, AS, and reference of origin are all reported
for the top two aligning references.

`{reference}-unique-raw.csv`: Contains data by species for reads that uniquely aligned to 


## Project Directory Structure

### `bin/`

A directory containing external scripts and/or binaries used for this pipeline. Currently only contains `assign-reads.py`, 
an external python script that extracts alignment scores for each intermediate bam produced in alignment.

*Note: this script relies on the python/3.12.3 cluster module having the polars library installed. Currently this is the case,
but ideally this would be changed to a conda environment or better, a docker container for robustness.*

### `data/`

A general purpose directory for storing data related to this project, whether it be supplemental files or a reference genome 
directory. 

If the `--univec_fasta` argument is not provided on pipeline execution, the univec fasta file will be downloaded to 
`data/univec/UniVec_Core.fa` and indexed with `bwa index` in that location further downstream.


