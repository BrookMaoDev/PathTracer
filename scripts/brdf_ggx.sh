#!/bin/sh
# GGX microfacet metal BRDF on the Cornell box, with IS and ES: the mirror sphere becomes gold
# (F0 = 1.0, 0.71, 0.29) with a configurable roughness.
#   ggx_gold_polished:  roughness 0.08 - sharp reflections and highlight
#   ggx_gold_burnished: roughness 0.3  - blurred reflections, broad highlight
. "$(dirname "$0")/common.sh"

build ggx_gold_polished -D__USE_IS -D__USE_ES -D__DEMO_GGX -DGGX_DEMO_ROUGHNESS=0.08
render ggx_gold_polished 512 7 1000

build ggx_gold_burnished -D__USE_IS -D__USE_ES -D__DEMO_GGX -DGGX_DEMO_ROUGHNESS=0.3
render ggx_gold_burnished 512 7 1000
