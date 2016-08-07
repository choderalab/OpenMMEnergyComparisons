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
pmdparm = pmd.load_file('step3_pbcsetup.psf')
pmdparm.positions = crd.positions
pmdparm.box = [80, 80, 80, 90, 90, 90]

omm_e = pmd.openmm.energy_decomposition_system(pmdparm, system)

print('OpenMM Energy is %s' % omm_e)
