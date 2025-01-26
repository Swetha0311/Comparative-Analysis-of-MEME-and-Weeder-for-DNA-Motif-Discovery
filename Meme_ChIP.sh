#!/bin/bash

# Meme on ChIP-seq data
meme /N/slate/sweyadav/motif/alignment/SRR392333_peak_sequences.fa -minw 6 -maxw 20 -dna -oc /N/slate/sweyadav/motif/SRR392333_meme

meme /N/slate/sweyadav/motif/alignment/SRR392334_peak_sequences.fa -minw 6 -maxw 20 -dna -oc /N/slate/sweyadav/motif/SRR392334_meme
