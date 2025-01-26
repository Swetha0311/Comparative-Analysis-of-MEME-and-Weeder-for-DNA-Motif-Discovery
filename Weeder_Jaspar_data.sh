#!/bin/bash

# Load necessary modules (if required)
module load intel/19.1.2

#Move to the directory
cd /N/slate/sweyadav/motif

# Paths
INPUT_DIR="/N/slate/sweyadav/motif/weeder_jas/jaspar_fixed"
WEEDER_BINARY="/N/slate/sweyadav/motif/weeder2"
TIME_LOG="$INPUT_DIR/time_memory_log.txt"

# Generate a list of input files
FILES=($(find "$INPUT_DIR" -type f -name "*.fa"))
FILE=${FILES[$SLURM_ARRAY_TASK_ID]} # Select the file based on the job array
index
BASENAME=$(basename "$FILE" .fa)

# Start time tracking
START_TIME=$(date +%s)

# Run Weeder with memory and time logging
/usr/bin/time -v "$WEEDER_BINARY" -f "$FILE" -O MM 2>
"$INPUT_DIR/${BASENAME}_time_mem.log"

# End time tracking
END_TIME=$(date +%s)
ELAPSED_TIME=$((END_TIME - START_TIME))

# Extract memory usage from /usr/bin/time output
MEMORY_USAGE=$(grep "Maximum resident set size"
"$INPUT_DIR/${BASENAME}_time_mem.log" | awk '{print $6}')

# Log the time and memory usage for each file
echo "File: $FILE, Time: ${ELAPSED_TIME}s, Memory: ${MEMORY_USAGE}KB" >>
"$TIME_LOG"