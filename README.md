# rhapsody-merge-microbial
Analysis code used for merging taxonomic classification results with BD Rhapsody outputs

## Background

Single cell sequencing data is typically used to characterize cells on the basis of their transcriptomic profiles.
However, sequencing data can also be used to characterize the microbial content of a sample.
While sequencing approaches are designed to only capture eukaryotic mRNA, microbial RNA/DNA may be incidentally captured and sequenced.

## Approach

After single-cell sequencing data has been generated, the raw sequencing reads can be processed
using existing tools for taxonomic classification, such as Kraken2.
Any such tool which provided a per-read classification can be used.
The code in this repository may be used to merge the taxonomic classification results with the BD Rhapsody outputs.
