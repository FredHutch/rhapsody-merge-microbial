#!/usr/bin/env python3

from collections import defaultdict
import pandas as pd
import gzip
import click


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


def map_tax_to_cell(tax_long: pd.DataFrame, fastq_1, fastq_2):

    # Build a table of taxonomic classification results
    output = defaultdict(lambda: defaultdict(int))

    # Iterate over both FASTQ files in parallel
    for (header1, seq1), (header2, seq2) in zip(
        parse_fastq(fastq_1),
        parse_fastq(fastq_2)
    ):
        # If there is no taxonomic classification result, skip it
        if tax_long.get(header1) is None:
            continue

        # Otherwise, get the cell barcode
        cell_barcode = get_cell_barcode(header1, seq1, header2, seq2)

        # Add the taxonomic classification result to the table
        output[cell_barcode][tax_long[header1]] += 1

    # Convert the table to a wide format
    tax_wide = pd.DataFrame(output).fillna(0)
    # report the number of reads classified to each taxon for each cell barcode
    print(f"Taxonomic classification results:\n{tax_wide.shape[0]:,} rows x {tax_wide.shape[1]:,} columns")
    return tax_wide


def get_cell_barcode(header1, seq1, header2, seq2):
    """Get the cell barcode for a given pair of reads."""
    return "BARCODE"  # FIXME


@click.command()
@click.argument('kraken_output', type=click.Path(exists=True))
@click.argument('barcodes_tsv', type=click.Path(exists=True))
@click.argument('bam', type=click.Path(exists=True))
def main(kraken_output, barcodes_tsv, bam):

    # Example usage:
    print(f"Kraken output file: {kraken_output}")
    print(f"Single cell table file: {barcodes_tsv}")
    print(f"BAM file: {bam}")

    # # Read in the table of taxonomic classification results
    # # Transform to a dict with read_id as key and tax_id as value
    # tax_long = pd.concat([
    #     (
    #         df
    #         .query("flag == 'C'")
    #         .query("tax_id != 1")
    #         .query("tax_id != 2")
    #         .query("tax_id != 9606")
    #         .query("tax_id != 28384")
    #         .query("tax_id != 131567")
    #     )
    #     for df in pd.read_csv(
    #         kraken_output,
    #         sep="\t",
    #         header=None,
    #         names=["flag", "read_id", "tax_id"],
    #         usecols=[0, 1, 2],
    #         iterator=True,
    #         chunksize=10000
    #     )
    # ]).set_index("read_id")["tax_id"].to_dict()
    # print(f"Taxonomic classification results for {len(tax_long):,} reads")
    # print("Top hits:")
    # print(pd.Series(tax_long).value_counts().head(10))

    # # Map each taxonomic classification result to the cell barcode
    # tax_wide = map_tax_to_cell(tax_long, fastq_1, fastq_2)

    # # Merge the table with the single cell table
    # # FIXME


if __name__ == '__main__':
    main()
