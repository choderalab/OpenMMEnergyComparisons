#!/usr/bin/env python

import parmed as pmd
from simtk import openmm as mm, unit as u
from simtk.openmm import app

psf = app.CharmmPsfFile('step3_pbcsetup.psf')
# Taken from output of CHARMM run
psf.setBox(6.3*u.nanometers, 6.3*u.nanometers, 6.3*u.nanometers)
crd = app.CharmmCrdFile('step3_pbcsetup.crd')
params = app.CharmmParameterSet('toppar/par_all36_na.prm',
                                'toppar/toppar_water_ions.str')

system = psf.createSystem(params, nonbondedMethod=app.PME,
        nonbondedCutoff=12*u.angstroms, switchDistance=10*u.angstroms)

for force in system.getForces():
    if isinstance(force, mm.CustomNonbondedForce):
        #print('CustomNonbondedForce: %s' % force.getUseSwitchingFunction())
        #print('LRC? %s' % force.getUseLongRangeCorrection())
        force.setUseLongRangeCorrection(False)
    elif isinstance(force, mm.NonbondedForce):
        #print('NonbondedForce: %s' % force.getUseSwitchingFunction())
        #print('LRC? %s' % force.getUseDispersionCorrection())
        force.setUseDispersionCorrection(False)
pmdparm = pmd.load_file('step3_pbcsetup.psf')
pmdparm.positions = crd.positions
pmdparm.box = [63, 63, 63, 90, 90, 90]

# CHARMM energy: From docker evaluation
# TODO: Pull these components in with a Python script.
"""
ENER ENR:  Eval#     ENERgy      Delta-E         GRMS
ENER INTERN:          BONDs       ANGLes       UREY-b    DIHEdrals    IMPRopers
ENER EXTERN:        VDWaals         ELEC       HBONds          ASP         USER
ENER IMAGES:        IMNBvdw       IMELec       IMHBnd       RXNField    EXTElec
ENER EWALD:          EWKSum       EWSElf       EWEXcl       EWQCor       EWUTil
 ----------       ---------    ---------    ---------    ---------    ---------
ENER>        0 -81613.68849      0.00000      5.11347
ENER INTERN>     3220.73435   2153.87430     86.79342    473.71332      0.91605
ENER EXTERN>    10734.99706 -81661.12324      0.00000      0.00000      0.00000
ENER IMAGES>      165.68305  -3570.40512      0.00000      0.00000      0.00000
ENER EWALD>      1706.56070-512279.28667 497353.85428      0.00000      0.00000
 ----------       ---------    ---------    ---------    ---------    ---------
"""
charmm_energy = dict()
charmm_energy['Bond'] = 3220.73435 
charmm_energy['Angle'] = 2153.87430    + 86.79342
charmm_energy['Dihedral'] = 473.71332 + 0.91605
charmm_energy['Nonbonded'] = 10734.99706 -81661.12324 + 165.68305  -3570.40512 + 1706.56070 -512279.28667 + 497353.85428
charmm_energy['Total'] = -81613.68849

total = 0.0
for key in ['Bond', 'Angle', 'Dihedral', 'Nonbonded']:
    total += charmm_energy[key]
print(charmm_energy['Total'], total)
omm_e = pmd.openmm.energy_decomposition_system(pmdparm, system)

# OpenMM energy
openmm_energy = dict()
openmm_energy['Bond'] = omm_e[0][1]
openmm_energy['Angle'] = omm_e[1][1] + omm_e[2][1]
openmm_energy['Dihedral'] = omm_e[3][1] + omm_e[4][1] + omm_e[5][1]
openmm_energy['Nonbonded'] = omm_e[6][1] + omm_e[7][1]
openmm_energy['Total'] = openmm_energy['Bond'] + openmm_energy['Angle'] + openmm_energy['Dihedral'] + openmm_energy['Nonbonded']

print('OpenMM Energy is %s' % omm_e)

# Now do the comparisons
print('Output in kJ/mol')
print('%-20s | %-15s | %-15s' % ('Component', 'CHARMM', 'OpenMM'))
print('-'*56)
total = 0
for name in ['Bond', 'Angle', 'Dihedral', 'Nonbonded']:
    print('%-20s | %15.6f | %15.6f' % (name, charmm_energy[name], openmm_energy[name]))
print('-'*56)
print('%-20s | %15.6f | %15.6f' % ('Total', charmm_energy['Total'] * 4.184, openmm_energy['Total'] * 4.184))
print('-'*56)
