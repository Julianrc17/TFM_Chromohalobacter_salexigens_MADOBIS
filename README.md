# TFM Chromohalobacter salexigens MADOBIS

## Project Description 

This repository contains the code and data used for my Master's Thesis, which focuses on a RNA-Seq of the halophilic bacterium Chromohalobacter salexigens. The primary goal of the study is to investigate the molecular mechanisms of salt tolerance in this organism. 
The project involves processing biological data, running sequence alignments, and performing statistical analyses to understand the underlying molecular pathways. Below you'll find detailed instructions on how to reproduce the analyses. 

I used the Linux environment on Windows through WSL (Windows Subsystem for Linux) to conduct RNA-Seq analysis. Starting with the raw data in .fastq format, I processed it step by step, ultimately generating count data in an Excel file. 

#### Analysis Data

- We have 14 samples with four conditions: Wt06, Wt25, RpoS06, and RpoS25.
- 06 corresponds to 0.6 M NaCl, considered low salinity, while 25 corresponds to 2.5 M NaCl, considered high salinity. All samples are at 37ºC.
- The samples obtained using the SOLID sequencing method are: Wt06A, Wt06B, Wt06C, Wt25A, Wt25B, and Wt25C. These files are in .bam format.
- The samples obtained using the Illumina sequencing method are: Wt06I, Wt25I, RpoS06A, RpoS06B, RpoS06C, RpoS25A, RpoS25B, and RpoS25C. These files are in .fastq format.


### RNA-Seq_RawData_Code.txt (Bash WSL)

This will be the first file encountered in the repository: a txt script that outlines, step by step, the code used to process the raw, unfiltered data. Based on a series of commands executed in WSL.

- cat command is used to join several files, in this case we have a series of .fastq files that are separate and to process and align them it is better to join them.
- fastqc quality control checks on the combined FASTQ files for read 1 and read 2, generating quality reports. All quality reports were checked to ensure they were of good quality and a cleaning was also done with trimmomatic but it is an optional step that in the end was not done in the final result.
- bwa mem aligns paired-end reads to the reference genome and outputs in SAM format.
- samtools view converts the SAM file to BAM format, and samtools sort sorts the BAM file.
- samtools flagstat provides summary statistics of the BAM file.
- samtools stats generates detailed alignment statistics.
