#!/bin/bash

# Load necessary modules
module load meme

# Paths
INPUT_DIR="/N/slate/sweyadav/motif/weeder_jas/jaspar_fixed"
OUTPUT_BASE_DIR="/N/slate/sweyadav/motif/meme_jas" # Base output directory
TIME_LOG="$OUTPUT_BASE_DIR/time_memory_log.txt"

# Ensure the base output directory exists
mkdir -p "$OUTPUT_BASE_DIR"

# Generate a list of input files
FILES=($(find "$INPUT_DIR" -type f -name "*.fa"))
FILE=${FILES[$SLURM_ARRAY_TASK_ID]} # Select the file based on the job array
index
BASENAME=$(basename "$FILE" .fa)

# Define the output directory for this file
OUTPUT_DIR="$OUTPUT_BASE_DIR/${BASENAME}_result"
mkdir -p "$OUTPUT_DIR" # Ensure the individual output directory exists

# Start time tracking
START_TIME=$(date +%s)

# Run MEME with input and output parameters
/usr/bin/time -v meme "$FILE" -oc "$OUTPUT_DIR" 2>
"$OUTPUT_DIR/${BASENAME}_time_mem.log"

# End time tracking
END_TIME=$(date +%s)
ELAPSED_TIME=$((END_TIME - START_TIME))

# Extract memory usage from /usr/bin/time output
MEMORY_USAGE=$(grep "Maximum resident set size"
"$OUTPUT_DIR/${BASENAME}_time_mem.log" | awk '{print $6}')

# Log the time and memory usage for each file
echo "File: $FILE, Time: ${ELAPSED_TIME}s, Memory: ${MEMORY_USAGE}KB" >>
"$TIME_LOG"