#!/bin/sh
# Cornell box, 512x512, 1000 samples, depth 7, no importance sampling, no explicit light sampling
. "$(dirname "$0")/common.sh"

build Cornell_noIS_noES
render Cornell_noIS_noES 512 7 1000
