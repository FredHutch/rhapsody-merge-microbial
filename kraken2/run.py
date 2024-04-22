#!/usr/bin/env python3

from collections import defaultdict
import pandas as pd
import gzip
import click
import pysam


def parse_fastq(fastq):
    with gzip.open(fastq, 'rt') as f:
        header = None
        for i, line in enumerate(f):
            if i % 4 == 0:
                assert line.startswith("@"), f"Unexpected header line: {line}"
                header = line[1:].strip().split(" ")[0]
            elif i % 4 == 1:
                yield header, line.strip()
            else:
                pass


def map_tax_to_cell(tax_long: pd.DataFrame, bam: str):

    # Make a lookup set with the read IDs which have been classified
    classified_reads = set(tax_long.index)

    print(f"Searching for cell barcodes for {len(classified_reads):,} reads in {bam}")

    # Make a dict with the cell barcode for each read ID
    cb = {}

    # Open the BAM file and iterate over every line
    with pysam.AlignmentFile(bam, "rb") as f:
        for read in f:
            # Skip reads which have not been classified
            if read.query_name not in classified_reads:
                continue

            # Extract the cell barcode from the read
            cell_barcode = get_cell_barcode(read)

            # Add the cell barcode to the taxonomic classification table
            cb[read.query_name, "cell_barcode"] = cell_barcode

    cb = pd.Series(cb)

    print(f"Found cell barcodes for {cb.shape[0]:,} reads")


def get_cell_barcode(read: pysam.AlignedSegment):
    """Get the cell barcode for a given pair of reads."""
    return read.get_tag("CB")


def get_tax_long(kraken_output):
    tax_long = pd.concat([
        (
            df
            .query("flag == 'C'")
            .query("tax_id != 1")
            .query("tax_id != 2")
            .query("tax_id != 9606")
            .query("tax_id != 28384")
            .query("tax_id != 131567")
        )
        for df in pd.read_csv(
            kraken_output,
            sep="\t",
            header=None,
            names=["flag", "read_id", "tax_id"],
            usecols=[0, 1, 2],
            iterator=True,
            chunksize=10000
        )
    ]).set_index("read_id")["tax_id"].to_dict()
    print(f"Taxonomic classification results for {len(tax_long):,} reads")
    print("Top hits:")
    print(pd.Series(tax_long).value_counts().head(10))

    return tax_long


@click.command()
@click.argument('kraken_output', type=click.Path(exists=True))
@click.argument('barcodes_tsv', type=click.Path(exists=True))
@click.argument('bam', type=click.Path(exists=True))
def main(kraken_output, barcodes_tsv, bam):

    # Example usage:
    print(f"Kraken output file: {kraken_output}")
    print(f"Single cell table file: {barcodes_tsv}")
    print(f"BAM file: {bam}")

    # Read in the table of taxonomic classification results
    # Transform to a dict with read_id as key and tax_id as value
    tax_long = get_tax_long(kraken_output)

    # Map each read ID to the cell barcode
    tax_wide = map_tax_to_cell(tax_long, bam)

    # # Merge the table with the single cell table
    # # FIXME


if __name__ == '__main__':
    main()
