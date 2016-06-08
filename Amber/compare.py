#!/usr/bin/env python

# For the OpenMM energies
import simtk.openmm as mm
import simtk.openmm.app as app
import simtk.unit as u

# For the Amber energies
import parmed as pmd
import sander

import argparse
import collections
import numpy as np

parser = argparse.ArgumentParser()

group = parser.add_argument_group('Input files')
group.add_argument('-p', '--prmtop', dest='prmtop', metavar='<FILE>',
                    required=True, help='Prmtop file to calculate energy')
group.add_argument('-c', '--inpcrd', dest='inpcrd', metavar='<FILE>',
                   required=True, help='Inpcrd file with coordinates (and box)')

group = parser.add_argument_group('Potential Energy Function')
group.add_argument('--pbc', dest='pbc', default=False, action='store_true',
                   help='Specify the use of periodic boundary conditions')
group.add_argument('-g', '--igb', dest='igb', default=0, type=int,
                   help='''GB model to use (see Amber manual). Default is
                   %(default)s''')
group.add_argument('--cutoff', dest='cutoff', type=float, default=8,
                   help='''Cutoff for nonbonded forces in Angstroms. Cutoff is
                   always infinite for non-periodic systems''')
group.add_argument('--platform', dest='platform', default='Reference',
                   help='OpenMM platform to use. Default is %(default)s')

args = parser.parse_args()

# Get the Amber forces and energies
if args.pbc:
    inp = sander.pme_input()
    inp.cut = args.cutoff
else:
    inp = sander.gas_input(args.igb)
    inp.cut = 1000

with sander.setup(args.prmtop, args.inpcrd, None, inp) as context:
    e, f = sander.energy_forces()
    f = np.array(f).reshape((context.natom, 3))

# Get the OpenMM forces and energies
parm = app.AmberPrmtopFile(args.prmtop)
inpcrd = app.AmberInpcrdFile(args.inpcrd)

gbmap = collections.defaultdict(lambda: None)
gbmap[1] = app.HCT
gbmap[2] = app.OBC1
gbmap[5] = app.OBC2
gbmap[7] = app.GBn
gbmap[8] = app.GBn2

if args.pbc:
    nonbondedMethod = app.PME
else:
    nonbondedMethod = app.NoCutoff

system = parm.createSystem(nonbondedMethod=nonbondedMethod,
                           nonbondedCutoff=args.cutoff*u.angstroms,
                           implicitSolvent=gbmap[args.igb])

pmdparm = pmd.load_file(args.prmtop, args.inpcrd)
omm_e = pmd.openmm.energy_decomposition_system(pmdparm, system,
                                               platform=args.platform)

# Generate the context and get the Force
con = mm.Context(system, mm.VerletIntegrator(1*u.femtoseconds),
                 mm.Platform.getPlatformByName(args.platform))
if args.pbc:
    con.setPeriodicBoxVectors(*inpcrd.box_vectors)
con.setPositions(inpcrd.positions)

omm_f = con.getState(getForces=True).getForces(asNumpy=True).value_in_unit(u.kilocalories_per_mole/u.angstroms)

# Now do the comparisons
print('%-20s | %-15s | %-15s' % ('Component', 'OpenMM', 'Amber'))
print('-'*56)
for name, oe in omm_e:
    if name == 'HarmonicBondForce':
        term, ae = 'Bond', e.bond
    elif name == 'HarmonicAngleForce':
        term, ae = 'Angle', e.angle
    elif name == 'PeriodicTorsionForce':
        term, ae = 'Dihedral', e.dihedral
    elif name == 'NonbondedForce':
        term, ae = 'Nonbonded', e.elec+e.vdw+e.elec_14+e.vdw_14
    elif name == 'CustomGBForce':
        term, ae = 'GB', e.egb
    else:
        continue
    print('%-20s | %15.6f | %15.6f' % (term, ae, oe))

# Now compare forces
proj = (f * omm_f).sum(axis=1) / (omm_f * omm_f).sum(axis=1)

print('Projection of Amber and OpenMM force:')
print('Average: %15.6f' % proj.mean())
print('Min:     %15.6f' % proj.min())
print('Max:     %15.6f' % proj.max())
