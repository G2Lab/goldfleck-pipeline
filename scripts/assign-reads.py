#!/usr/bin/env python3

import argparse
import polars as pl
import pysam

from functools import reduce
from pathlib import Path


"""
Assign reads from a set of BAM files (same reads, aligned to different references)
to a reference genome. All bam files are parsed and all output data can be traced 
back to said bam file by its stem (filename with extension stripped).

By: Akash Pandit
Last edited: July 31st 2025

Arguments:
-d | --bamdir <filepath>: path to directory containing bam files. 
-o | --outdir <filepath>: path to output directory, default="./output-assign-reads"

Outputs:

outdir/top-alignments.csv
    csv containing the top 2 AS / NM / MAPQ scores for all reads in which both
    top 2 species pass QC, along with the reference bamfile name.

outdir/{species}-unique-raw.csv
    csv containing unique read IDs and AS / MAPQ / NM for reads that uniquely
    mapped + passed QC for this species. One file is generated per bamfile in bamdir.

outdir/{species}-unique-top-align.csv
    csv containing top alignments (as seen in top-alignments.csv) for reads that 
    uniquely mapped + passed QC for this species and 2nd place. One file is generated
    per bamfile in bamdir.
"""


def load_bamfile_df(bampath: str | Path) -> pl.DataFrame:
    """ 
    Extracts read ID (qname + _R1/2), mapq, NM, and AS values for a bamfile (only
    for primary/unmapped read entries) and dumps them in a polars DataFrame.
    """

    if not isinstance(bampath, Path):
        bampath = Path(str(bampath))
    species = bampath.stem
    
    with pysam.AlignmentFile(str(bampath), "rb") as bamfile:
        records = [
            (
                f"{read.query_name}_{'R1' if read.is_read1 else 'R2'}",  # qname + _R1/2, unique ID
                read.mapping_quality,                                    # mapq
                int(read.get_tag("NM")) if read.has_tag("NM") else -1,   # num mismatches, opt sam tag
                int(read.get_tag("AS")) if read.has_tag("AS") else 0,    # align score, opt sam tag, bwa mem defaults 0
            )
            for read in bamfile.fetch(until_eof=True)
            if not (read.is_secondary or read.is_supplementary)
        ]

    return (
        pl.DataFrame(
            data=records,
            orient="row",
            schema={
                "unique_read_id": pl.String,
                f"{species}_mapq": pl.Int16,
                f"{species}_nm": pl.Int16,
                f"{species}_as": pl.Int16,
            }
        )
        .with_columns(pl.col(f"{species}_nm").fill_null(-1))
    )


def filter_for_mapq_nm(species: str) -> pl.Expr:
    """ 
    creates and returns the filter (mapq in [20, 30) and NM < 6) or (mapq 30+ and
    NM < 4) on columns for the given species.
    """
    return ( 
        pl.col(species+"_mapq").is_between(20, 30, closed="left")   # mapq in [20, 30)
      & pl.col(species+"_nm").lt(6)                                 # and nm < 6
    ) | (                                                           # or
        pl.col(species+"_mapq").ge(30)                              # mapq 30+
      & pl.col(species+"_nm").lt(4)                                 # and nm < 4
    )


def get_top_2(
    row: tuple,
    as_species: list[tuple[str, int]],
    nm_indices: dict[str, int],
    mapq_indices: dict[str, int]
) -> tuple:
    sample = row[0]

    # Build list of (species, score) using precomputed AS positions
    as_scores = [(species, row[idx]) for species, idx in as_species]

    if len(as_scores) < 2:
        return (sample,) + (None,) * 6  # Pad with Nones if fewer than 2 scores

    # Sort by alignment score and take top 2
    top2 = sorted(as_scores, key=lambda x: x[1], reverse=True)[:2]
    (top_fp, top_as), (second_fp, second_as) = top2

    # Lookup NM / MAPQ values using precomputed column indices
    top_nm    = row[nm_indices[top_fp]]
    second_nm = row[nm_indices[second_fp]]
    top_mapq    = row[mapq_indices[top_fp]]
    second_mapq = row[mapq_indices[top_fp]]

    return sample, top_as, second_as, top_nm, second_nm, top_mapq, second_mapq, top_fp, second_fp


if __name__ == "__main__":
    """
    1. handle CLI arguments
    """


    parser = argparse.ArgumentParser()
    parser.add_argument("-d", "--bamdir", type=Path, action="store", required=True)
    parser.add_argument("-o", "--outdir", type=Path, action="store", required=False, default=Path("output-assign-reads"))
    args = parser.parse_args()
    bamdir: Path = args.bamdir
    outdir: Path = args.outdir

    # argument validation
    if not bamdir.is_dir():
        print("Error: path", str(bamdir), "is not an existing directory.")
        exit(1)

    outdir.mkdir(exist_ok=True)


    """
    2. get merged dataframe w/ qc & as fields
    """


    species_paths = list(bamdir.glob("*.bam"))
    species_list  = [path.stem for path in species_paths]

    source_df = load_bamfile_df(species_paths[0])
    print("[Status] loaded bamfile into dataframe for", species_paths[0].name)

    for species_path in species_paths[1:]:
        species_df = load_bamfile_df(species_path)
        print("[Status] loaded bamfile into dataframe for", species_path.name)
        source_df = source_df.join(other=species_df, on="unique_read_id")
        del species_df  # every million reads adds ~51Mb ram used, wise to clean up large objects


    """
    3. dump rows that uniquely meet qc conditions for a species into csvs
    """


    passing_ids: dict[str, set] = {}

    for species in species_list:
        qc_expr = filter_for_mapq_nm(species)
        species_df = source_df.filter(qc_expr)
        passing_ids[species] = set(species_df.get_column("unique_read_id"))

    unique_ids: dict[str, set] = {}

    for species, idset in passing_ids.items():
        complement_set = set().union(*[st for sp, st in passing_ids.items() if sp != species])
        unique_ids[species] = idset.difference(complement_set)

    for species, uidset in unique_ids.items():
        cols = ["unique_read_id", species + "_as", species + "_mapq", species + "_nm"]
        (
            source_df
            .select(cols)
            .filter(pl.col("unique_read_id").is_in(uidset))
            .write_csv(outdir/(species + "-unique-raw.csv"))
        )
        print("[Status] dumped unique read csv for", species)


    """
    4. remove rows which don't meet qc conditions for any species
    """


    # as long as qc passes for any species, the row is kept
    qc_conditions = [filter_for_mapq_nm(species) for species in species_list]
    qc_conditions = reduce(lambda cond1, cond2: cond1 | cond2, qc_conditions)

    # set up some constants 
    columns = source_df.columns
    as_species  =  [(col[:-3], i) for i, col in enumerate(columns) if col.endswith("_as")]
    nm_indices   = { col[:-3]: i  for i, col in enumerate(columns) if col.endswith("_nm")}
    mapq_indices = { col[:-5]: i  for i, col in enumerate(columns) if col.endswith("_mapq")}

    top_as,   second_as   = pl.col("top_as"),   pl.col("second_as")
    top_mapq, second_mapq = pl.col("top_mapq"), pl.col("second_mapq")
    top_nm,   second_nm   = pl.col("top_nm"),   pl.col("second_nm")

    print("[Status] Extracting top 2 alignments for read set...")
    filtered_df = (
        source_df
        .filter(qc_conditions)
        .map_rows(lambda row: get_top_2(row, as_species, nm_indices, mapq_indices))
        .rename(mapping={
            "column_0": "unique_read_id",
            "column_1": "top_as",
            "column_2": "second_as",
            "column_3": "top_nm",
            "column_4": "second_nm",
            "column_5": "top_mapq",
            "column_6": "second_mapq",
            "column_7": "top_species",
            "column_8": "second_species"
        })
        .filter(
            (
                ( top_mapq.is_between(20, 29) & top_nm.lt(6) )
                |
                ( top_mapq.ge(30) & top_nm.lt(4) )
            )
            &
            (
                ( second_mapq.is_between(20, 29) & second_nm.lt(6) )
                |
                ( second_mapq.ge(30) & second_nm.lt(4) )
            )
        )
        .filter(top_as.ne(second_as))
    )
    filtered_df.write_csv(outdir/"top-alignments.csv")


    """
    5. dump top alignments for unique reads by species
    """

    for species, uidset in unique_ids.items():
        (
            filtered_df
            .filter(pl.col("unique_read_id").is_in(uidset))
            .write_csv(outdir/(species + "-unique-top-alignments.csv"))
        )
        print("[Status] dumped unique read top alignments csv for", species)
