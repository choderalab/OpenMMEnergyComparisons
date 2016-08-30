#!/usr/bin/env python

import parmed as pmd
from simtk import openmm as mm, unit as u
from simtk.openmm import app

psf = app.CharmmPsfFile('step3_pbcsetup.psf')
# Taken from output of CHARMM run
psf.setBox(8*u.nanometers, 8*u.nanometers, 8*u.nanometers)
crd = app.CharmmCrdFile('step3_pbcsetup.crd')
params = app.CharmmParameterSet('toppar/par_all36_prot.prm',
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
        force.setPMEParameters(1.0/0.34, 90, 90, 90)
pmdparm = pmd.load_file('step3_pbcsetup.psf')
pmdparm.positions = crd.positions
pmdparm.box = [80, 80, 80, 90, 90, 90]

# Print PME parameters
integrator = mm.VerletIntegrator(1.0 * u.femtoseconds)
context = mm.Context(system, integrator)
context.setPositions(crd.positions)
for force in system.getForces():
    if isinstance(force, mm.NonbondedForce):
        break
integrator.step(1)
print(force.getPMEParametersInContext(context))
del context, integrator

# CHARMM energy: From docker evaluation
# TODO: Pull these components in with a Python script.
"""
ENER ENR:  Eval#     ENERgy      Delta-E         GRMS
ENER INTERN:          BONDs       ANGLes       UREY-b    DIHEdrals    IMPRopers
ENER CROSS:           CMAPs        PMF1D        PMF2D        PRIMO
ENER EXTERN:        VDWaals         ELEC       HBONds          ASP         USER
ENER IMAGES:        IMNBvdw       IMELec       IMHBnd       RXNField    EXTElec
ENER EWALD:          EWKSum       EWSElf       EWEXcl       EWQCor       EWUTil
 ----------       ---------    ---------    ---------    ---------    ---------
ENER>        0-163043.99835      0.00000      5.02279
ENER INTERN>     6337.99813   4236.12181     54.30685   1726.66813     21.86301
ENER CROSS>       -21.48984      0.00000      0.00000      0.00000
ENER EXTERN>    20161.20647-164737.82886      0.00000      0.00000      0.00000
ENER IMAGES>      243.39096  -5318.48694      0.00000      0.00000      0.00000
ENER EWALD>       4130.5989-1021718.0599  991839.7129       0.0000       0.0000
"""
charmm_energy = dict()
charmm_energy['Bond'] = 6337.99813
charmm_energy['Angle'] = 4236.12181 + 54.30685
charmm_energy['Dihedral'] = 1726.66813 + 21.86301 + -21.48984
charmm_energy['Nonbonded'] = 20161.20647 + -164737.82886 + 243.39096  -5318.48694 + 4130.5989 -1021718.0599 + 991839.7129
charmm_energy['Total'] = -163043.99835

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
    print('%-20s | %15.2f | %15.2f' % (name, charmm_energy[name] * 4.184, openmm_energy[name] * 4.184))
print('-'*56)
print('%-20s | %15.2f | %15.2f' % ('Total', charmm_energy['Total'] * 4.184, openmm_energy['Total'] * 4.184))
print('-'*56)
