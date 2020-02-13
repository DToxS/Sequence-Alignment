# Sequence Alignment

This analysis uses the *STAR* sequence aligner to map all the sequence reads contained in the *FASTQ* data files obtained from high-throughput mRNA-sequencing assay to the UCSC human genome reference library (hg38), such that those uniquely aligned reads can be assigned to corresponding reference genes in the <u>Feature Counts</u> analysis later.

## Inputs

The inputs of this analysis include:

- The gzipped *FASTQ* data files in the `Seqs` directory.
- The UCSC human genome reference library in the `References` directory.

## Procedure

The procedure of this analysis includes the following steps:

1. Build a set of index files for the human genome reference library used by *STAR*, by running `Programs/Sequence-Alignment/Build-STAR-Genome-Index.GEO.sh`.

2. Align a set of mRNA-sequencing *FASTQ* data files to the human genome reference library, by running `Programs/Sequence-Alignment/Align-Conv-RNAseq-Reads.GEO.sh`.


## Outputs

The outputs of this analysis include:

- The *STAR* genome index data directory `STAR-Index-[Version]` in the `References` directory.
- A set of sequence alignment data files (**.bam*) in the `Aligns` directory.

