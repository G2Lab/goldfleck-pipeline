# Multi-genome Alignment Pipeline

## !! NOTE !!: currently outdated, following changes in place:
- `wgsim` used to simulate 1 paired-end read set with a user-defined blend of reference sequence : proportion (percent) pairs
- analysis moved from basic analysis of flagstat output to more comprehensive analysis of map quality, mismatch number, and alignment score from bamfiles
- WIP: porting to single nextflow script for ease of use
- TBD: rewrite this document for said nextflow script

A pipeline designed to:
- Simulate a set of paired-end reads with varying mutation rates from a set of reference genomes using `wgsim`
- Align simulated reads to each reference genome and generate `bam` files for each pair of reads, filtering out known contaminants (e.g. PCR primers) using NCBI's UniVec database
- Analyze bamfile statistics and demonstrate that bamfile from a given species `A` can be uniquely mapped back to species `A`, and not the other species

## Scripts

### `scripts/filter-and-align.sh`

Aligns provided FASTQ files to every genome (FASTA file) provided, first discarding known contaminants in [NCBI's UniVec database](https://www.ncbi.nlm.nih.gov/tools/vecscreen/univec/). FASTQs are aligned to each genome in parallel with generated SLURM scripts. Generated scripts can be found with their outputs in the output directory.

Usage: `sbatch filter-and-align.sh -1 <read1.fastq.gz> -2 <read2.fastq.gz> -d <refdirs> [-u <univec>] [-o <outdir>]`

| Argument | Description |
|-------------------|-------------------------|
| `-1 <file.fastq>` | paired-end fastq file 1 |
| `-2 <file.fastq>` | paired-end fastq file 2 |
| `-d <directory>` | reference sequences (structured as directory/genome/genome.fna) |
| `-u <univec path>` | (optional) path to a univec reference file. Downloaded from NCBI FTP to current directory if not found |
| `-o <outdir>` | (optional) path to output directory, `align-output` if unspecified |
| `-h` | show help message and exit, ignores other arguments |

#### Output file structure:

```
align-out/
├── bamfiles/
│   ├── refseq-filename.bam
│   ... (for each refseq)
│
├── fastqs/
│   ├── paired1.fq  (univec filtered fastqs)
│   └── paired2.fq
│
└── slurm/
    ├── slurm-[job ID].out
    ├── slurm-align-species.sh
    ... (for each refseq)
```


### `scripts/gen-wgsim-reads.py`

Generates simulated FASTQs from all reference genomes provided. For each reference, FASTQs are produced with mutation rates (`-r` argument) 0, 0.001, 0.0025, 0.005, 0.01, 0.025, and 0.05. Each species is run in parallel with a generated slurm script, with each script running wgsim with all mutation rates. Other wgsim arguments are maintained as follows:
- `-e` base error rate: 0.020
- `-d` outer distance: 100 bp
- `-s` standard deviation: 30 bp
- `-N` number of read pairs: 1,000,000
- `-1` length of first read: 150
- `-2` length of second read: 150

Usage: `python3 gen-wgsim-reads.py [-h] -d <path> [-o <path>] [-s <seed>]`

| Argument | Description |
|-------------------|-------------------------|
| `-d/--ref-dir <path>` | reference sequences (structured as directory/genome/genome.fna) |
| `-o/--out-dir <path>` | (optional) path to output directory, `wgsim-output` if unspecified |
| `-s/--seed <seed>` | (optional) random seed for wgsim generation, only needed for exact reproducability |
| `-h/--help` | show help message and exit, ignores other arguments |

#### Output file structure:

```
wgsim-out/
├── fastqs/
│   ├── refseq/
│   │   ├── logs/
│   │   │   ├── slurm-refseq.out
│   │   │   ├── wgsim-refseqmut[r].log  (note: logs include data on inserted mutations)
│   │   │   ... (for each mutation rate)
│   │   │
│   │   ├── refseq-mut[r]-1.fq
│   │   ├── refseq-mut[r]-2.fq
│   │   ... (for each mutation rate)
│   │
│   ... (for each refseq)
│
├── refseq-inputs/ (symbolic link that points to actual -d location)
│
└── scripts/
    ├── wgsim-slurm-species.sh
    ... (for each refseq)
```

### `scripts/run-alignment.sh`

**Not a main part of the pipeline.** 

Batch runs `filter-and-align.sh` for every mutation rate only using lion generated FASTQs. Takes no arguments.

Usage: `bash run-alignment.sh` 

