#!/bin/sh
# Texture, normal, and alpha mapping (ported from A3) on the Cornell box, with IS and ES:
#   left sphere:  stone colour texture + normal map
#   right sphere: gold paint colour texture + lattice alpha map (rays and shadow rays pass through the holes)
. "$(dirname "$0")/common.sh"

build texture_normal_alpha_mapping -D__USE_IS -D__USE_ES -D__DEMO_MAPPING
render texture_normal_alpha_mapping 512 7 1000
