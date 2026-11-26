#!/bin/sh
# Shared setup for the render scripts (sourced, not run directly).
#   build NAME [FLAGS...]               compiles the path tracer with the given -D flags into bin/NAME
#   render NAME SIZE DEPTH SAMPLES      renders renderings/NAME.ppm at SIZE x SIZE with that binary
set -e
cd "$(dirname "$0")/.."	# Repo root
mkdir -p bin renderings
export OMP_WAIT_POLICY=passive	# Idle threads sleep at the per-pass barrier instead of spinning (~10% faster here)

build() {
  name=$1
  shift
  g++ -O3 -fopenmp src/svdDynamic.c src/PathTracer.c src/utils_path.c -lm "$@" -o "bin/$name"
}

render() {
  start=$(date +%s)
  "./bin/$1" "$2" "$3" "$4" "renderings/$1.ppm"
  echo "$1: rendered renderings/$1.ppm in $(( $(date +%s) - start ))s"
}
