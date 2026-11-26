#!/bin/sh
# Cornell box, 512x512, 1000 samples, depth 7, importance sampling and explicit light sampling
. "$(dirname "$0")/common.sh"

build Cornell_IS_ES -D__USE_IS -D__USE_ES
render Cornell_IS_ES 512 7 1000
