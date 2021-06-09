# Calling PacBio Methylation 
Small pipeline to call methylation in Streptococcus suis using standard PacBio tools

## Goal

This pipeline was constructed to call m6A and m4C methylation in *Streptococcus suis* mutants. These mutants should differ in their methylation status, which we wanted to investigate using PacBio. For an introduction on identification of methylation using PacBio data, see [this Nature Methods paper](https://dx.doi.org/10.1038%2Fnmeth.1459).

:warning: The pipeline was uploaded to GitHub to be able to retrace our steps and as a rough guide for future analyses. This pipeline is not intended as a tool which others can easily apply to their own data :warning: 

## Requirements

This pipeline has the following requirements:
- [Snakemake](https://snakemake.readthedocs.io/) and [conda](https://docs.conda.io/en/latest/) are available on your system
- ipdSummary and MultiMotifMaker are installed on your system. This is not possible through conda. See https://www.pacb.com/support/software-downloads/ for ipdSummary (part of smrttools) and https://github.com/bioinfomaticsCSU/MultiMotifMaker for MultiMotifMaker. Update paths to executables as appropriate.
- PacBio read data is available in fasta format in the directory `pacbio_fa` and paired-end Illumina data is available in gzipped fastq format in the directory `illumina_fastq`
- Enough computational power. The analysis was originally run on a 96 GB RAM computing node with 16 CPUs. This worked alright.

Unicycler and pbmm2 will be installed through conda, integrated in Snakemake.

## Methods

This pipeline combines several tools to end up with motifs associated with particular methylation. The tools and their functions are:

| Tool | Purpose | Outputs |
|----|----|----|
| Unicycler | Hybrid assembly of complete genome | Assembly fasta (`asembly.fasta`) and some quality control files |
| pbmm2 | Mapping of PacBio reads to assembly | Bam file ("native PacBio" format) |
| ipdSummary | Identification of m6A, m4C and unknown modifications from mapped reads | Modifications file in GFF3 format and an exhaustive overview of kinetics per nt in csv format |
| MultiMotifMaker | Identification of motifs associated with modifications | Csv file containing motifs |

The final csv file containing motifs is the main output.
