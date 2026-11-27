#!/bin/sh
# Depth of field (thin lens camera) on the Cornell box, with IS and ES.
#   dof_focus:          circular aperture, focused on the glass sphere; mirror sphere and back wall blur
#   dof_bokeh_hexagon:  wide 6-blade aperture focused in front of the spheres; highlights blur into hexagons
. "$(dirname "$0")/common.sh"

build dof_focus -D__USE_IS -D__USE_ES -D__USE_DOF -DDOF_APERTURE=1 -DDOF_FOCUS=19.5
render dof_focus 512 7 1000

build dof_bokeh_hexagon -D__USE_IS -D__USE_ES -D__USE_DOF -DDOF_APERTURE=3 -DDOF_FOCUS=15.5 -DDOF_BLADES=6
render dof_bokeh_hexagon 512 7 1000
