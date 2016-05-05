#!/bin/sh
# Make sure we don't go past any errors
set -e

echo "Downloading the PDB files 4m6j and 2koc"
python << EOF
import parmed as pmd
pmd.download_PDB('4m6j', saveto='4m6j.pdb')
pmd.download_PDB('2koc', saveto='2koc.pdb')
EOF

# First strip down the PDB file
echo "Fixing the PDB file (adding H, removing water, adjusting residue names)"
pdb4amber -i 4m6j.pdb -o 4m6j_fixed.pdb --reduce --dry
pdb4amber -i 2koc.pdb -o 2koc_fixed.pdb --dry --model=1 # NMR structure already has Hs

echo "Deleting the cofactor in DHFR and 5' phosphate in hairpin"
# Delete the NDP residue
python << EOF
import parmed as pmd
pmd.load_file('4m6j_fixed.pdb')['!:NDP'].save('4m6j_fixed.pdb', overwrite=True)
pmd.load_file('2koc_fixed.pdb')[4:].save('2koc_fixed.pdb', overwrite=True)
EOF

# Now create the prmtop
echo "Creating the DHFR topology file"
tleap -f - << EOF
source leaprc.protein.ff14SB
source leaprc.water.tip4pew
sys = loadPDB 4m6j_fixed.pdb
set default pbradii mbondi3
saveAmberParm sys dhfr.dry.parm7 dhfr.dry.rst7
solvateOct sys TIP4PEWBOX 15.0
addIonsRand sys NA 20 CL 20
saveAmberParm sys dhfr.pbc.parm7 dhfr.pbc.rst7
quit
EOF

echo "Creating the hairpin RNA topology file"
tleap -f - << EOF
source leaprc.RNA.OL3
source leaprc.water.tip4pew
sys = loadPDB 2koc_fixed.pdb
set default pbradii mbondi3
saveAmberParm sys 2koc.dry.parm7 2koc.dry.rst7
solvateOct sys TIP4PEWBOX 15.0
addIonsRand sys NA 20 CL 20
saveAmberParm sys 2koc.pbc.parm7 2koc.pbc.rst7
quit
EOF
