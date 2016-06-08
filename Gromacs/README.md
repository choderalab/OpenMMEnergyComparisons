GROMACS files
=============

What we do here is just copy the systems from CHARMM and Amber to do more kinds
of comparisons.

CHARMM Systems
--------------

CHARMM-GUI already comes with GROMACS inputs, so those files
should be copied here directly.

AMBER Systems
-------------

The AMBER files are converted using ParmEd.  The ``convert.sh`` file will run
``convert_amber.py`` twice on the two solvated systems in the Amber directory to
generate the GROMACS topology and GRO files here.
