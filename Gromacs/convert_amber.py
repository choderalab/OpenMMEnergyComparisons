#!/usr/bin/env python

import parmed as pmd
import sys

if len(sys.argv) != 5:
    sys.exit('Usage: %s <prmtop> <inpcrd> <gromacs_top> <gromacs_gro>' %
             sys.argv[0])

prmtop, inpcrd, gmxtop, gmxgro = sys.argv[1:]

parm = pmd.load_file(prmtop, inpcrd)
parm.save(gmxtop, format='gromacs')
parm.save(gmxgro, format='gro')
