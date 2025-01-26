# Comparative-Analysis-of-MEME-and-Weeder-for-DNA-Motif-Discovery

# Overview

This project compares the performance of two DNA motif discovery tools, MEME and Weeder, using ChIP-Seq-derived data and JASPAR motifs. The study evaluates their efficiency, accuracy, and ability to detect biologically significant motifs.

# Project Goals

Analyze MEME and Weeder for motif discovery in Mus musculus sequences

Compare performance based on computational efficiency and accuracy

Evaluate motif similarity with known JASPAR motifs

# Datasets Used

JASPAR Database – 254 FASTA files containing 3,845,662 sequences

ChIP-Seq Data (PRJNA142963) – GSM851274 (RenLab-H3K27ac-kidney)

SRR392333: 479.5M bases, 13,319,058 reads

SRR392334: 450.9M bases, 12,524,950 reads

# Methodology

1. Preprocessing

JASPAR: Standardized headers for compatibility

ChIP-Seq: Quality control (FASTQC), adapter trimming (TrimGalore), genome alignment (STAR), peak calling (HOMER)

2. Motif Discovery

MEME: Expectation-Maximization (EM) algorithm for motif identification

Weeder: Exhaustive search for shorter conserved motifs

3. Evaluation

Accuracy: Comparison with JASPAR motifs using TOMTOM

Computational Performance: Time and memory usage benchmarking

# Key Findings

MEME efficiently detects longer motifs but struggles with alternative binding motifs

Weeder identifies shorter conserved motifs with higher accuracy but has high computational demands

MEME is significantly faster, while Weeder is more biologically accurate for ChIP-Seq data

# Dependencies

MEME Suite

Weeder

STAR Aligner

HOMER

TrimGalore

FASTQC

R (ggplot2 for visualization)
