# Sequence Alignment

This step of analysis uses the *STAR* sequence aligner to map all the sequence reads contained in the *FASTQ* data files obtained from high-throughput mRNA-sequencing assay to the UCSC human genome reference library (*hg38*), such that those uniquely aligned reads can be assigned to corresponding reference genes in the [Feature Counts](https://github.com/DToxS/Feature-Counts) analysis later.

**Note:** the `[Type]` tag included in the file and directory names below refers to either `Conv` for conventional sequence data or `DGE` for 3'-DGE sequence data.

## Inputs

The inputs of this analysis include:

- The gzipped *FASTQ* data files in the `Seqs` sub-directory of each dataset directory under the `[Type]-GEO-Depot` top directory.
- The UCSC human genome reference library in the `References` directory under the `[Type]-GEO-Depot` top directory.

## Procedure

The procedure of this analysis includes the following steps:

1. Set the variable `DATASET_DIR` in `Programs/Sequence-Alignment/Build-[Type]-STAR-Genome-Index.GEO.sh` to the absolute path of the `[Type]-GEO-Depot` top directory, and launch the program to build a set of index files for the human genome reference library used by *STAR*.
2. Set the variable `DATASET_DIR` in `Programs/Sequence-Alignment/Align-[Type]-RNAseq-Reads.GEO.sh` to the absolute path of the `[Type]-GEO-Depot` top directory, and launch the program to align corresponding `[Type]` of mRNA-sequencing *FASTQ* data files to the human genome reference library.


## Outputs

The outputs of this analysis include:

- The *STAR* genome index data directory `STAR-Index-[Version]-[Type]` in the `References` directory.
- A set of sequence alignment data files (`*.bam`) in the `Aligns` sub-directory of each dataset directory under the corresponding `[Type]-GEO-Depot` top directory.

