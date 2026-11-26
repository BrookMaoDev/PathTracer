#!/bin/sh
# Cornell box, 512x512, 1000 samples, depth 7, importance sampling, no explicit light sampling
. "$(dirname "$0")/common.sh"

build Cornell_IS_noES -D__USE_IS
render Cornell_IS_noES 512 7 1000
