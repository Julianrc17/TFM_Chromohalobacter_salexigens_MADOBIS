# Complete Code for RNA-Seq

### Analysis Data
We have 14 samples with four conditions: Wt06, Wt25, RpoS06, and RpoS25.
06 corresponds to 0.6 M NaCl, considered low salinity, while 25 corresponds to 2.5 M NaCl, considered high salinity. All samples are at 37ºC.
The samples obtained using the SOLID sequencing method are: Wt06A, Wt06B, Wt06C, Wt25A, Wt25B, and Wt25C.
The samples obtained using the Illumina sequencing method are: Wt06I, Wt25I, RpoS06A, RpoS06B, RpoS06C, RpoS25A, RpoS25B, and RpoS25C.

For the analysis, we will also need the .gtf and .fna files, which can be downloaded from the following NCBI page where the latest updates of the Chromohalobacter salexigens genome are available: https://www.ncbi.nlm.nih.gov/datasets/taxonomy/158080/.
The RefSeq is: GCF_000055785.1

## ILLUMINA SAMPLES
bwa index GCF_000055785.1_ASM5578v1_genomic.fna
