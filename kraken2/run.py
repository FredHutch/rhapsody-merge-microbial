#!/usr/bin/env python3

from collections import defaultdict
import pandas as pd
import gzip
import click
import pysam
import os
from functools import lru_cache
from collections import defaultdict
import logging


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


def map_tax_to_cell(tax_long: dict, bam: str):

    print(f"Searching for cell barcodes for {len(tax_long):,} reads in {bam}")

    # Make a dict with the cell barcode for each read ID
    cb = {}

    # Open the BAM file and iterate over every line
    with pysam.AlignmentFile(bam, "rb") as f:
        for read in f:
            # Skip reads which have not been classified
            if read.query_name not in tax_long:
                continue

            # Extract the cell barcode from the read
            cell_barcode = get_cell_barcode(read)

            # Add the cell barcode to the taxonomic classification table
            cb[read.query_name] = cell_barcode

            if len(cb) >= 10:
                print("FIXME -- stopping the BAM parsing early")
                break

    print(f"Found cell barcodes for {len(cb):,} reads")

    # Make a wide DataFrame with CB as the index and tax_id as the columns
    tax_wide = pd.DataFrame(dict(
        cb=cb,
        tax=tax_long
    ))

    tax_wide = tax_wide.pivot_table(
        index="cb",
        columns="tax",
        aggfunc=len
    )

    return tax_wide


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
            chunksize=10000,
            nrows=10000
        )
    ]).set_index("read_id")["tax_id"].to_dict()
    print(f"Taxonomic classification results for {len(tax_long):,} reads")
    print("Top hits:")
    print(pd.Series(tax_long).value_counts().head(10))

    return tax_long


def make_tax_summary(tax_wide: pd.DataFrame):

    tax = NCBITaxonomy('ncbi_taxonomy')
    return (
        pd.DataFrame([
            tax.summary(str(tax_id))
            for tax_id in tax_wide.columns
        ])
        .set_index("tax_id")
        .assign(
            total_classified_reads=tax_wide.sum()
        )
    )


# Read in the taxonomy
class NCBITaxonomy():
    def __init__(self, folder):
        self.tax = defaultdict(dict)
        # Read in the file of taxid information
        names_fp = os.path.join(folder, 'names.dmp')
        assert os.path.exists(names_fp)
        with open(names_fp) as f:
            for line in f:
                line = line.strip('\n').split('\t')
                taxid = line[0]
                val = line[2]
                key = line[6]
                assert key != '', line
                self.tax[taxid][key] = val
        # Read in the file linking taxids to rank and parents
        nodes_fp = os.path.join(folder, 'nodes.dmp')
        assert os.path.exists(nodes_fp)
        with open(nodes_fp) as f:
            for line in f:
                line = line.strip('\n').split('\t')
                child = line[0]
                parent = line[2]
                child_rank = line[4]
                assert parent in self.tax, line
                self.tax[child]['rank'] = child_rank
                if parent != child:
                    self.tax[child]['parent'] = parent

        # Read in the "merged" taxids
        merged_fp = os.path.join(folder, 'merged.dmp')
        assert os.path.exists(merged_fp)
        with open(merged_fp) as f:
            for line in f:
                line = line.strip('\n').split('\t')
                old_taxid = line[0]
                new_taxid = line[2]
                assert new_taxid in self.tax
                # Link the old taxid to the new taxid
                self.tax[old_taxid] = self.tax[new_taxid]

    def summary(self, taxid):
        return dict(
            tax_id=taxid,
            name=self.name(taxid),
            rank=self.rank(taxid),
            lineage=";".join([
                self.name(taxid)
                for taxid in self.path_to_root(taxid)[::-1]
            ])
        )

    def info(self, taxid):
        assert taxid in self.tax, "Tax ID not found: {}".format(taxid)
        return self.tax[taxid]

    def name(self, taxid):
        assert taxid in self.tax, "Tax ID not found: {}".format(taxid)
        return self.tax[taxid]['scientific name']

    def rank(self, taxid):
        assert taxid in self.tax, "Tax ID not found: {}".format(taxid)
        return self.tax[taxid]['rank']

    def path_to_root(self, taxid):
        visited = [taxid]
        while 'parent' in self.tax[taxid]:
            taxid = self.tax[taxid]['parent']
            if 'ancestors' in self.tax[taxid]:
                visited.extend(self.tax[taxid]['ancestors'])
                return visited
            assert taxid not in visited, visited
            visited.append(taxid)
        return visited

    def is_below(self, taxid, group_taxid):
        # Determine whether a taxid is part of a particular group
        assert taxid in self.tax, "Tax ID not found: {}".format(taxid)
        if 'ancestors' not in self.tax[taxid]:
            self.tax[taxid]['ancestors'] = self.path_to_root(taxid)

        return group_taxid in self.tax[taxid]['ancestors']

    @lru_cache(maxsize=None)
    def lca(self, taxid1, taxid2):
        # Return the lowest common ancestor of both taxid1 and taxid2
        # Get the list of ancestors for both taxids (starting at the root)
        anc1 = self.path_to_root(taxid1)[::-1]
        anc2 = self.path_to_root(taxid2)[::-1]
        # Both lists end at the root
        if anc1[0] != anc2[0]:
            logging.info("{} and {} not rooted on the same taxonomy, returning None".format(taxid1, taxid2))
            return None
        # Set the initial LCA as the root
        lca = anc1[0]

        # Walk down the list of ancestors until they no longer agree
        for a1, a2 in zip(anc1, anc2):
            # If the ancestors are in common, update the LCA
            if a1 == a2:
                lca = a1
            else:
                break

        return lca

    @lru_cache(maxsize=None)
    def anc_at_rank(self, taxid, rank):
        """Return the ancestor of this taxid at a specific rank."""
        # The result is potentially None (if `taxid` is above `rank`)
        if self.rank(taxid) == rank:
            return taxid
        for t in self.path_to_root(taxid):
            if self.rank(t) == rank:
                return t
        return None

    def domain(self, taxid):
        """Determine the broad organismal domain that a taxid belongs to."""
        ancestors = self.path_to_root(taxid)
        if "2" in ancestors:
            d = "Bacteria"
        elif "10239" in ancestors:
            d = "Viruses"
        elif "2157" in ancestors:
            d = "Archaea"
        elif "4751" in ancestors:
            d = "Fungi"
        elif "33208" in ancestors:
            d = "Metazoa"
        elif "33090" in ancestors:
            d = "Green plants"
        elif "2759" in ancestors:
            d = "Other Eukaryotes"
        else:
            d = "Non-microbial"
        return d
   
    def maximum_weighted_path(self, list_of_taxids):
        """Get the maximum weighted path for a list of taxids."""
        taxid_weights = defaultdict(int)
        for taxid in list_of_taxids:
            for anc in self.path_to_root(taxid):
                if anc in list_of_taxids:
                    taxid_weights[taxid] += 1
        max_count = max(taxid_weights.values())
        mwp = [k for k, v in taxid_weights.items() if v == max_count]
        
        # If there is more than one MWP, return the LCA of them
        while len(mwp) > 1:
            mwp[1] = self.lca(mwp[0], mwp[1])
            mwp = mwp[1:]
            
        return mwp[0]


@click.command()
@click.argument('kraken_output', type=click.Path(exists=True))
@click.argument('bam', type=click.Path(exists=True))
@click.argument('output_prefix', type=str)
def main(kraken_output, bam, output_prefix):

    # Example usage:
    print(f"Kraken output file: {kraken_output}")
    print(f"BAM file: {bam}")

    # Read in the table of taxonomic classification results
    # Transform to a dict with read_id as key and tax_id as value
    tax_long = get_tax_long(kraken_output)

    # Map each read ID to the cell barcode
    tax_wide = map_tax_to_cell(tax_long, bam)

    # Write out the wide table
    tax_wide.to_csv(f"{output_prefix}.counts_per_cell.csv")

    # Make a summary for each taxon which includes the name, lineage, and total reads
    tax_summary = make_tax_summary(tax_wide)


if __name__ == '__main__':
    main()
