#!/bin/bash

# Load the SRA Toolkit for downloading SRA data
module load sra-toolkit

# Download sequencing data from SRA using prefetch for multiple SRA IDs
prefetch SRR392333 SRR392334

# Load FastQC module for quality control checks
module load fastqc

# Convert SRA files to FASTQ format using fasterq-dump
fasterq-dump SRR392333 SRR392334

# Run FastQC quality control on each FASTQ file
fastqc SRR392333 SRR392334

# Load Trim Galore module for adapter trimming and quality control
module load trimgalore

# Perform trimming and generate quality control reports
trim_galore --adapter GATCGGAAGAGCACACGTCTGAACTCCAGTCACATC --fastqc --stringency 3 SRR392333.fastq

# Load the STAR module for genome indexing
module load star/2.7.11a

# Genome indexing
STAR --runMode genomeGenerate --genomeDir /N/slate/sweyadav/motif/index --genomeFastaFiles /N/slate/sweyadav/motif/GRCh38.p14.genome.fa --sjdbGTFfile /N/slate/sweyadav/motif/gencode.v44.chr_patch_hapl_scaff.annotation.gtf --sjdbOverhang 100

# Genome alignment
for FILENAME in SRR392333 SRR392334; do
STAR \
--outSAMattributes All \
--outSAMtype BAM SortedByCoordinate \
--quantMode GeneCounts \
--readFilesCommand cat \
--runThreadN 2 \
--sjdbGTFfile /N/slate/sweyadav/motif/gencode.vM36.chr_patch_hapl_scaff.annotation.gtf \
--outReadsUnmapped Fastx \
--outMultimapperOrder Random \
--outWigType wiggle \
--genomeDir /N/slate/sweyadav/motif/index \
--readFilesIn /N/slate/sweyadav/motif/${FILENAME}_trimmed.fq \
--outFileNamePrefix /N/slate/sweyadav/motif/alignment/${FILENAME}_
done

# Get peak files
module load bedtools/2.31.0
bedtools getfasta -fi /N/slate/sweyadav/motif/GRCm39.genome.fa -bed SRR392334_peaks.bed -fo SRR392334_peak_sequences.fa
