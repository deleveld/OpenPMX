#!/bin/bash

# Exit immediately if a command exits with a non-zero status
set -e

echo "Starting cleanup of specified file extensions..."

# Find and remove files matching the given extensions recursively
find . -type f \( \
    -name "*.out" -o \
    -name "*.ext" -o \
    -name "*.cov" -o \
    -name "*.yhat" -o \
    -name "*.phi" -o \
    -name "*.iter" -o \
    -name "*.gr.c" -o \
    -name "*.covcheck" -o \
    -name "*.o" \
\) -exec rm -v {} +

echo "Cleanup complete!"
