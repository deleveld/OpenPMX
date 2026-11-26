#!/bin/bash

# Exit immediately if a command exits with a non-zero status
set -e

echo "Starting cleanup of specified file extensions..."

# Find and remove files matching the given extensions recursively
find . -type f \( \
    -name "nonmem" -o \
    -name "*.out" -o \
    -name "*.ext" -o \
    -name "*.cov" -o \
    -name "*.yhat" -o \
    -name "*.phi" -o \
    -name "*.png" -o \
    -name "*.iter" -o \
    -name "*.icov" -o \
    -name "*.gr.c" -o \
    -name "*.gr.graphics" -o \
    -name "*.gr.pdf" -o \
    -name "*.gr.pdf" -o \
    -name "*.bin" -o \
    -name "*.cov.limits" -o \
    -name "*.profie.iter" -o \
    -name "*.vpc.points" -o \
    -name "*.o" \
\) -exec rm -v {} +

echo "Cleanup complete!"
