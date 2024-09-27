# TFM Chromohalobacter salexigens MADOBIS

## Project Description 

This repository contains the code and data used for my Master's Thesis, which focuses on a RNA-Seq of the halophilic bacterium Chromohalobacter salexigens. The primary goal of the study is to investigate the molecular mechanisms of salt tolerance in this organism. 
The project involves processing biological data, running sequence alignments, and performing statistical analyses to understand the underlying molecular pathways. Below you'll find detailed instructions on how to reproduce the analyses. 

I used the Linux environment on Windows through WSL (Windows Subsystem for Linux) to conduct RNA-Seq analysis. Starting with the raw data in .fastq format, I processed it step by step, ultimately generating count data in an Excel file where the rest of the analysis was carried out. 

#### Analysis Data

- We have 14 samples with four conditions: Wt06, Wt25, RpoS06, and RpoS25.
- 06 corresponds to 0.6 M NaCl, considered low salinity, while 25 corresponds to 2.5 M NaCl, considered high salinity. All samples are at 37ºC.
- The samples obtained using the SOLID sequencing method are: Wt06A, Wt06B, Wt06C, Wt25A, Wt25B, and Wt25C. These files are in .bam format.
- The samples obtained using the Illumina sequencing method are: Wt06I, Wt25I, RpoS06A, RpoS06B, RpoS06C, RpoS25A, RpoS25B, and RpoS25C. These files are in .fastq format.


### RNA-Seq_RawData_Code.txt (Bash WSL)

This will be the first file encountered in the repository: a txt script that outlines, step by step, the code used to process the raw, unfiltered data. Based on a series of commands executed in WSL.

- **cat** command is used to join several files, in this case we have a series of .fastq files that are separate and to process and align them it is better to join them.
- **fastqc** quality control checks on the combined FASTQ files for read 1 and read 2, generating quality reports. All quality reports were checked to ensure they were of good quality and a cleaning was also done with trimmomatic but it is an optional step that in the final result was not done.
- **bwa mem** aligns paired-end reads to the reference genome and outputs in SAM format.
- **samtools view** converts the SAM file to BAM format, and samtools sort sorts the BAM file.
- **samtools flagstat** provides summary statistics of the BAM file.
- **samtools stats** generates detailed alignment statistics.

**Program required**: **Bwa mem** (Burrows-Wheeler Alignment Tool), **FastQC**, **Samtools**. 

### R Analysis 

On the other hand, we have the analysis done in R. Which deals with everything from the processing and normalization of our count data obtained with FeatureCounts to the generation of FoldChange tables and relevance graphs in our study.

**Data Loading & Preparation:**
Raw count data from FeatureCounts is loaded and preprocessed, including renaming columns and creating a DESeqDataSet object for analysis.

**Data Normalization:**
Geometric means are used to normalize the counts, and size factors are estimated to account for sequencing depth differences across samples.

**Differential Expression Analysis:**
DESeq2 package from BioConductor is used to compare conditions, computing log2 Fold Changes. A custom function transforms log2 Fold Changes to regular Fold Changes, and results are filtered by significance (p-value and FoldChange thresholds). Filtered results are written to Excel files for each comparison.

**Visualization:**
Bar plots show the number of over-expressed and under-expressed genes at different FoldChange thresholds. Volcano plots visualize significant gene expression changes, highlighting relevant genes for further study.

**The rest of the analysis will be carried out using Excel to process the data**
