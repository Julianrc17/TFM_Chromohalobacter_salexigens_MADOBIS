#!/bin/bash

# Complete Code for RNA-Seq

# Analysis Data

# We have 14 samples with four conditions: Wt06, Wt25, RpoS06, and RpoS25.
# 06 corresponds to 0.6 M NaCl, considered low salinity, while 25 corresponds to 2.5 M NaCl, considered high salinity. All samples are at 37ºC.
# The samples obtained using the SOLID sequencing method are: Wt06A, Wt06B, Wt06C, Wt25A, Wt25B, and Wt25C.
# The samples obtained using the Illumina sequencing method are: Wt06I, Wt25I, RpoS06A, RpoS06B, RpoS06C, RpoS25A, RpoS25B, and RpoS25C.

# For the analysis, we will also need the .gtf and .fna files, which can be downloaded from the following NCBI page where the latest updates of the Chromohalobacter salexigens genome are available: 
# https://www.ncbi.nlm.nih.gov/datasets/taxonomy/158080/.
# The RefSeq is: GCF_000055785.1

# General Process of analysis described in README file.

# ILLUMINA SAMPLES

# Index of reference genome

bwa index GCF_000055785.1_ASM5578v1_genomic.fna

# Wt06I

cat Wt06_ATGTCA_L008_R1_001.fastq Wt06_ATGTCA_L008_R1_002.fastq Wt06_ATGTCA_L008_R1_003.fastq Wt06_ATGTCA_L008_R1_004.fastq > combinedWt06_R1.fastq
cat Wt06_ATGTCA_L008_R2_001.fastq Wt06_ATGTCA_L008_R2_002.fastq Wt06_ATGTCA_L008_R2_003.fastq Wt06_ATGTCA_L008_R2_004.fastq > combinedWt06_R2.fastq
../../FastQC./fastqc/ combinedWt06_R1.fastq combinedWt06_R2.fastq
bwa mem ../../Genoma_referencia_gtf/GCF_000055785.1_ASM5578v1_genomic.fna combinedWt06_R1.fastq combinedWt06_R2.fastq > aligned_reads_Wt06.sam
samtools view -Sb aligned_reads_Wt06.sam | samtools sort -o aligned_reads_Wt06_sorted.bam
samtools flagstat aligned_reads_Wt06_sorted.bam > aligned_stats_Wt06.txt
samtools stats aligned_reads_sorted.bam > alignment_stats_Wt06_full.txt

# Wt25I

cat Wt25_CCGTCC_L008_R1_001.fastq Wt25_CCGTCC_L008_R1_002.fastq Wt25_CCGTCC_L008_R1_003.fastq Wt25_CCGTCC_L008_R1_004.fastq Wt25_CCGTCC_L008_R1_005.fastq > combinedWt25_R1.fastq
cat Wt25_CCGTCC_L008_R2_001.fastq Wt25_CCGTCC_L008_R2_002.fastq Wt25_CCGTCC_L008_R2_003.fastq Wt25_CCGTCC_L008_R2_004.fastq Wt25_CCGTCC_L008_R2_005.fastq > combinedWt25_R2.fastq
../../FastQC./fastqc/ combinedWt25_R1.fastq combinedWt25_R2.fastq
bwa mem ../../Genoma_referencia_gtf/GCF_000055785.1_ASM5578v1_genomic.fna combinedWt25_R1.fastq combinedWt25_R2.fastq > aligned_reads_Wt25.sam
samtools view -Sb aligned_reads_Wt25.sam | samtools sort -o aligned_reads_Wt25_sorted.bam
samtools flagstat aligned_reads_Wt25_sorted.bam > aligned_stats_Wt25.txt
samtools stats aligned_reads_Wt25_sorted.bam > alignment_stats_Wt25_full.txt

# RpoS06A

cat RpoS06A_ATCACG_L008_R1_001.fastq RpoS06A_ATCACG_L008_R1_002.fastq RpoS06A_ATCACG_L008_R1_003.fastq > combinedRpos06A_R1.fastq
cat RpoS06A_ATCACG_L008_R2_001.fastq RpoS06A_ATCACG_L008_R2_002.fastq RpoS06A_ATCACG_L008_R2_003.fastq > combinedRpos06A_R2.fastq
../../FastQC./fastqc/ combinedRpoS06A_R1.fastq combinedRpoS06A_R2.fastq
bwa mem ../../Genoma_referencia_gtf/GCF_000055785.1_ASM5578v1_genomic.fna combinedRpoS06A_R1.fastq combinedRpoS06A_R2.fastq > aligned_reads_RpoS06A.sam
samtools view -Sb aligned_reads_RpoS06A.sam | samtools sort -o aligned_reads_RpoS06A_sorted.bam
samtools flagstat aligned_reads_RpoS06A_sorted.bam > aligned_stats_RpoS06A.txt
samtools stats aligned_reads_RpoS06A_sorted.bam > alignment_stats_RpoS06A_full.txt

# RpoS06B

cat RpoS06B_CGATGT_L008_R1_001.fastq RpoS06B_CGATGT_L008_R1_002.fastq RpoS06B_CGATGT_L008_R1_003.fastq RpoS06B_CGATGT_L008_R1_004.fastq RpoS06B_CGATGT_L008_R1_005.fastq > combinedRpoS06B_R1.fastq
cat RpoS06B_CGATGT_L008_R2_001.fastq RpoS06B_CGATGT_L008_R2_002.fastq RpoS06B_CGATGT_L008_R2_003.fastq RpoS06B_CGATGT_L008_R2_004.fastq RpoS06B_CGATGT_L008_R2_005.fastq > combinedRpoS06B_R2.fastq
../../FastQC./fastqc/ combinedRpoS06B_R1.fastq combinedRpoS06B_R2.fastq
bwa mem ../../Genoma_referencia_gtf/GCF_000055785.1_ASM5578v1_genomic.fna combinedRpoS06B_R1.fastq combinedRpoS06B_R2.fastq > aligned_reads_RpoS06B.sam
samtools view -Sb aligned_reads_RpoS06B.sam | samtools sort -o aligned_reads_RpoS06B_sorted.bam
samtools flagstat aligned_reads_RpoS06B_sorted.bam > aligned_stats_RpoS06B.txt
samtools stats aligned_reads_RpoS06B_sorted.bam > alignment_stats_RpoS06B_full.txt

# RpoS06C

cat RpoS06C_TTAGGC_L008_R1_001.fastq RpoS06C_TTAGGC_L008_R1_002.fastq RpoS06C_TTAGGC_L008_R1_003.fastq RpoS06C_TTAGGC_L008_R1_004.fastq > combinedRpoS06C_R1.fastq
cat RpoS06C_TTAGGC_L008_R2_001.fastq RpoS06C_TTAGGC_L008_R2_002.fastq RpoS06C_TTAGGC_L008_R2_003.fastq RpoS06C_TTAGGC_L008_R2_004.fastq > combinedRpoS06C_R2.fastq
../../FastQC./fastqc/ combinedRpoS06C_R1.fastq combinedRpoS06C_R2.fastq
bwa mem ../../Genoma_referencia_gtf/GCF_000055785.1_ASM5578v1_genomic.fna combinedRpoS06C_R1.fastq combinedRpoS06C_R2.fastq > aligned_reads_RpoS06C.sam
samtools view -Sb aligned_reads_RpoS06C.sam | samtools sort -o aligned_reads_RpoS06C_sorted.bam
samtools flagstat aligned_reads_RpoS06C_sorted.bam > aligned_stats_RpoS06C.txt
samtools stats aligned_reads_RpoS06C_sorted.bam > alignment_stats_RpoS06C_full.txt

# RpoS25A

cat RpoS25A_TGACCA_L008_R1_001.fastq RpoS25A_TGACCA_L008_R1_002.fastq RpoS25A_TGACCA_L008_R1_003.fastq RpoS25A_TGACCA_L008_R1_004.fastq RpoS25A_TGACCA_L008_R1_005.fastq > combinedRpoS25A_R1.fastq
cat RpoS25A_TGACCA_L008_R2_001.fastq RpoS25A_TGACCA_L008_R2_002.fastq RpoS25A_TGACCA_L008_R2_003.fastq RpoS25A_TGACCA_L008_R2_004.fastq RpoS25A_TGACCA_L008_R2_005.fastq > combinedRpoS25A_R2.fastq
../../FastQC./fastqc/ combinedRpoS25A_R1.fastq combinedRpoS25A_R2.fastq
bwa mem ../../Genoma_referencia_gtf/GCF_000055785.1_ASM5578v1_genomic.fna combinedRpoS25A_R1.fastq combinedRpoS25A_R2.fastq > aligned_reads_RpoS25A.sam
samtools view -Sb aligned_reads_RpoS25A.sam | samtools sort -o aligned_reads_RpoS25A_sorted.bam
samtools flagstat aligned_reads_RpoS25A_sorted.bam > aligned_stats_RpoS25A.txt
samtools stats aligned_reads_RpoS25A_sorted.bam > alignment_stats_RpoS25A_full.txt

# RpoS25B

cat RpoS25B_ACAGTG_L008_R1_001.fastq RpoS25B_ACAGTG_L008_R1_002.fastq RpoS25B_ACAGTG_L008_R1_003.fastq RpoS25B_ACAGTG_L008_R1_004.fastq > combinedRpoS25B_R1.fastq
cat RpoS25B_ACAGTG_L008_R2_001.fastq RpoS25B_ACAGTG_L008_R2_002.fastq RpoS25B_ACAGTG_L008_R2_003.fastq RpoS25B_ACAGTG_L008_R2_004.fastq > combinedRpoS25B_R2.fastq
../../FastQC./fastqc/ combinedRpoS25B_R1.fastq combinedRpoS25B_R2.fastq
bwa mem ../../Genoma_referencia_gtf/GCF_000055785.1_ASM5578v1_genomic.fna combinedRpoS25B_R1.fastq combinedRpoS25B_R2.fastq > aligned_reads_RpoS25B.sam
samtools view -Sb aligned_reads_RpoS25B.sam | samtools sort -o aligned_reads_RpoS25B_sorted.bam
samtools flagstat aligned_reads_RpoS25B_sorted.bam > aligned_stats_RpoS25B.txt
samtools stats aligned_reads_RpoS25B_sorted.bam > alignment_stats_RpoS25B_full.txt
 
# RpoS25C

cat RpoS25C_GCCAAT_L008_R1_001.fastq RpoS25C_GCCAAT_L008_R1_002.fastq RpoS25C_GCCAAT_L008_R1_003.fastq RpoS25C_GCCAAT_L008_R1_004.fastq > combinedRpoS25_R1.fastq
cat RpoS25C_GCCAAT_L008_R2_001.fastq RpoS25C_GCCAAT_L008_R2_002.fastq RpoS25C_GCCAAT_L008_R2_003.fastq RpoS25C_GCCAAT_L008_R2_004.fastq > combinedRpoS25C_R2.fastq
../../FastQC./fastqc/ combinedRpoS25C_R1.fastq combinedRpoS25C_R2.fastq
bwa mem ../../Genoma_referencia_gtf/GCF_000055785.1_ASM5578v1_genomic.fna combinedRpoS25C_R1.fastq combinedRpoS25C_R2.fastq > aligned_reads_RpoS25C.sam
samtools view -Sb aligned_reads_RpoS25C.sam | samtools sort -o aligned_reads_RpoS25C_sorted.bam
samtools flagstat aligned_reads_RpoS25C_sorted.bam > aligned_stats_RpoS25C.txt
samtools stats aligned_reads_RpoS25C_sorted.bam > alignment_stats_RpoS25C_full.txt

# SOLID SAMPLES

# These samples were already .bam. There is no need to pre-process them. Only statistical data is taken.

# Wt06A

samtools flagstat mergeBam_06A_Dr.bam > report_Wt06Abam.txt
samtools stats mergeBam_06A_Dr.bam > alignment_stats_06A_full.txt

# Wt06B
samtools flagstat mergeBam_06B_Dr.bam > report_Wt06Bbam.txt
samtools stats mergeBam_06B_Dr.bam > alignment_stats_06B_full.txt

# Wt06C
samtools flagstat mergeBam_06C_Dr.bam > report_W06Cbam.txt
samtools stats mergeBam_06C_Dr.bam > alignment_stats_06C_full.txt

# Wt25A
samtools flagstat mergeBam_25A_Dr.bam > report_Wt25Abam.txt
samtools stats mergeBam_25A_Dr.bam > alignment_stats_25A_full.txt

# Wt25B
samtools flagstat mergeBam_25B_Dr.bam > report_Wt25Bbam.txt
samtools stats mergeBam_25B_Dr.bam > alignment_stats_25B_full.txt

# Wt25C
samtools flagstat mergeBam_25C_Dr.bam > report_Wt25Cbam.txt
samtools stats mergeBam_25C_Dr.bam > alignment_stats_25C_full.txt

# FeatureCounts for All Samples

featureCounts -a ../../Genoma_referencia_gtf/genomic.gtf -o counts_totales.txt -t gene -g gene_id -p -T 12 ../../../bam\ solid/mergeBam_06A_Dr.bam ../../../bam\ solid/mergeBam_06B_Dr.bam ../../../bam\ solid/mergeBam_06C_Dr.bam ../../../bam\ solid/mergeBam_25A_Dr.bam ../../../bam\ solid/mergeBam_25B_Dr.bam ../../../bam\ solid/mergeBam_25C_Dr.bam bwa/aligned_reads_Wt25_sorted.bam ../Sample_Wt06/bwa/aligned_reads_Wt06_sorted.bam ../../RpoS/Sample_RpoS06A/bwa/aligned_reads_RpoS06A_sorted.bam ../../RpoS/Sample_RpoS06B/bwa/aligned_reads_RpoS06B_sorted.bam ../../RpoS/Sample_RpoS06C/bwa/aligned_reads_RpoS06C_sorted.bam ../../RpoS/Sample_RpoS25A/bwa/aligned_reads_RpoS25A_sorted.bam ../../RpoS/Sample_RpoS25B/bwa/aligned_reads_RpoS25B_sorted.bam ../../RpoS/Sample_RpoS25C/aligned_reads_RpoS25C_sorted.bam



