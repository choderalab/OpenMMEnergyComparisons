#!/bin/sh

# Convert all of the files we need to

python convert_amber.py ../Amber/2koc.pbc.parm7 ../Amber/2koc.pbc.rst7 2koc.pbc.top 2koc.pbc.gro
python convert_amber.py ../Amber/dhfr.pbc.parm7 ../Amber/dhfr.pbc.rst7 dhfr.pbc.top dhfr.pbc.gro
