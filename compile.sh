#!/bin/sh
# Run from the repo root. Sources live in src/, the binary is written to the repo root.
cd "$(dirname "$0")"

# Use this command to compile single-threaded for developing and debugging your code
#g++ -O3 -g src/svdDynamic.c src/PathTracer.c src/utils_path.c -lm -o PathTracer

# Use this command to compile multi-threaded so you can quickly render complex images
# once your code is working
g++ -O3 -g -fopenmp src/svdDynamic.c src/PathTracer.c src/utils_path.c -lm -o PathTracer
