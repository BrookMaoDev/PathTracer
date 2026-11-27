#!/bin/sh
# Spectral dispersion on the Cornell box, with IS and ES: the glass sphere uses Cauchy B = 0.15
# (exaggerated; real crown glass is ~0.0042), so its caustic and edges split into colours.
# Extra samples since each path carries a single wavelength (colour noise).
. "$(dirname "$0")/common.sh"

build dispersion -D__USE_IS -D__USE_ES -D__USE_DISPERSION -DDISPERSION_B=0.15
render dispersion 512 7 2000
