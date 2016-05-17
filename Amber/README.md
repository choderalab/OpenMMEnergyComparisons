Preparing the Amber input files
===============================

You need to have AmberTools 16 (released April 30, 2016) installed on your
computer for this to work. It uses ParmEd to download and process the PDB files
(in addition to pdb4amber), and automates the entire procedure of preparing the
prmtop and inpcrd files.

To generate the input files, use the command

```
./prepare.sh
```

This will generate the following files

- dhfr.gas.parm7
- dhfr.gas.rst7
- dhfr.pbc.parm7
- dhfr.pbc.rst7
- 2koc.gas.parm7
- 2koc.gas.rst7
- 2koc.pbc.parm7
- 2koc.pbc.rst7

These will be used to run the GB (``*.gas.*``) and PME (``*.pbc.*``)
comparisons, respectively.

Comparing Energies
==================

The Python script ``compare.py`` will evaluate energies using Amber and OpenMM
and put the results in a table.  Use the ``--help`` flag for help.
