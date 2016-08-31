CHARMM Input File Generation
============================

For the protein
---------------

To generate the input files, follow these steps:

1. Go to [CHARMM GUI's MD setup utility](http://www.charmm-gui.org/?doc=input/mdsetup)
2. Enter ``4mj6`` in the RCSB PDB download box and click ``Next Step``
3. Stick with the default selection (Protein, but no Hetero or Water) and click ``Next Step``
4. Select ``Terminal group patching`` and make sure the selected termini are
   ``NTER`` and ``CTER``. Then click ``Next Step``
5. Select ``Fit Waterbox Size to Protein Size`` and enter ``15.0`` in the
   ``Enter Edge Distance`` box. Keep remaining defaults (0.15 M KCl ions with
   Monte Carlo placement, which should add 43 K+ and 43 Cl- ions). Then click
   ``Next Step`` (you may have to re-enter ion information if requested)
6. Stick with the default of automatic grid generation and click ``Next Step``
7. Keep the default options in the next screen, generating input files for all
   available programs. Click ``Next Step``.
8. Finally, click the ``download .tgz`` button to download the generated files.

From the charmm-gui.tgz file, extract the ``toppar/`` directory, all of the
``step3_*`` files, ``namd/`` directory (if you plan on using NAMD to do energy
comparisons), and the ``openmm/`` directory.

For the RNA
-----------

We will follow basically the same steps as we did for the protein.  Naturally
there will be some changes, and only these are highlighted below.

2. Enter ``2koc`` instead of ``4mj6``
3. The default selection here will be RNA, and there will be no Hetero or Water.
4. Keep the selected ``Preserve hydrogen coordinates`` option (since this is an
   NMR structure with hydrogen positions). The termini should be 5TER and 3TER
   rather than NTER and CTER).

From the charmm-gui.tgz file, extract the ``toppar/`` directory, all of the
``step3_*`` files, ``namd/`` directory (if you plan on using NAMD to do energy
comparisons), and the ``openmm/`` directory.

CHARMM docker image
===================

Prerequisites
-------------
* [Docker toolbox](https://www.docker.com/products/docker-toolbox)
* [CHARMM lite nonprofit/academic version](http://charmm.chemistry.harvard.edu/charmm_lite.php) (free) - downloaded as `charmm.tar.gz`

Building the docker image
-------------------------
After starting the Docker daemon, run
```
docker build -t omnia/charmm-lite:c40b1 .
```

Running CHARMM
--------------
The CHARMM executable is `/charmm/c40b1_gnu/exec/gnu/charmm`

To manually start the docker image (for testing purposes):
```
docker run -i -t omnia/charmm-lite:c40b1 /bin/bash

Running comparison
------------------
```
To run the comparison from this directory:
```
python energy.py dhfr
python energy.py 2koc
```
