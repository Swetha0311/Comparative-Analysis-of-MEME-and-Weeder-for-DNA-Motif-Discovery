#!/bin/bash

# Load necessary modules (if required)
module load intel/19.1.2

# Paths
TIME_LOG="time_memory_log.txt"

# Start time tracking
START_TIME=$(date +%s)

#Run the code
./weeder2 -f /N/slate/sweyadav/motif/weeder/SRR392333_peak_sequences.fa -O MM

# End time tracking
END_TIME=$(date +%s)

# Calculate elapsed time
ELAPSED_TIME=$((END_TIME - START_TIME))

# Record memory usage (using SLURM sacct)
MEMORY_USAGE=$(sacct -j ${SLURM_JOB_ID} --format=MaxRSS --noheader | awk
'{print $1}')

# Log the time and memory usage
echo "Job: $SLURM_JOB_ID, Time: ${ELAPSED_TIME}s, Memory:
${MEMORY_USAGE}KB" >> $TIME_LOG

# Optional: Print out time and memory usage to the standard output (in addition to
logging it)
echo "Time: ${ELAPSED_TIME}s, Memory: ${MEMORY_USAGE}KB"