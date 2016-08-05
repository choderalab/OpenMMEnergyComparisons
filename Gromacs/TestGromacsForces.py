from __future__ import print_function
from simtk.openmm import *
from simtk.openmm.app import *
from simtk.unit import *
import os
import subprocess
import re

def printUsage():
    print("TestGromacsForces.py requires at least one argument of the form filename=(filename)")
    print("filename should be the name (without extension) of the top/gro/mdp files to load")
    print("It also accepts additional arguments of the form argumentName=value:")
    print("platform:  A platform (Reference, CPU, OpenCL, or CUDA); default is platform=Reference")
    print("precision: The precision to use (single, mixed, or double); default is single")
    print()
    print("Example:   python TestGromacsForces.py filename=dhfr.pbc platform=cuda")

# Set default values

filename = None
platformName = 'Reference'
precision = 'single'

# Parse the argument list

for arg in sys.argv[1:]:
    argSplit = arg.split("=")
    argName = argSplit[0].lower();
    argValue = argSplit[1]
    if argName=='filename':
        filename = argValue
    elif argName=='platform':
        platformName = argValue
    elif argName=='precision':
        precision = argValue

# Ensure that the argument filename is provided

if filename is None:
    print("Error: Argument list must contain a filename")
    printUsage()
    sys.exit(1) 

print("filename =", filename)
print("platformName =", platformName)
print("precision =", precision)

if platformName == 'OpenCL':
    properties = {'OpenCLPrecision':precision}
elif platformName == 'CUDA':
    properties = {'CudaPrecision':precision}
else:
    properties = dict()

# Compute forces with Gromacs

gromacsBinDir = '/usr/local/gromacs/bin/'
gromacsTopDir = '/usr/local/gromacs/share/gromacs/top'
subprocess.call([gromacsBinDir+'grompp', '-f', filename+'.mdp', '-p', filename+'.top', '-c', filename+'.gro'])
subprocess.call([gromacsBinDir+'mdrun', '-nt', '1'])
process = subprocess.Popen([gromacsBinDir+'gmxdump', '-f', 'traj.trr'], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
(output, errors) = process.communicate()
expr = re.compile('f\[.*?\]=\{(.*?),(.*?),(.*?)\}')
gromacsForces = []
for match in re.findall(expr, output):
    gromacsForces.append(Vec3(*[float(x) for x in match]))
terms = {}
terms['bonds'] = {'Bond'}
terms['angles'] = {'Angle'}
terms['torsions'] = {'Proper Dih.', 'Improper Dih.'}
terms['nonbonded'] = {'LJ-14', 'Coulomb-14', 'LJ (SR)', 'Coulomb (SR)', 'Coul. recip.', 'Disper. corr.'}
termEnergies = dict((term, 0.0) for term in terms)
process = subprocess.Popen([gromacsBinDir+'gmxdump', '-e', 'ener.edr'], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
(output, errors) = process.communicate()
for line in output.split('\n'):
    line = line.strip()
    for term in terms:
        for prefix in terms[term]:
            if line.startswith(prefix):
                termEnergies[term] += float(line[len(prefix):].split()[0])

# Compute forces with OpenMM

gro = GromacsGroFile(filename+'.gro')
top = GromacsTopFile(filename+'.top', periodicBoxVectors=gro.getPeriodicBoxVectors(), includeDir=gromacsTopDir)
mdpOptions = {}
for line in open(filename+'.mdp'):
    if '=' in line:
        fields = [x.strip().lower() for x in line.split('=')]
    mdpOptions[fields[0]] = fields[1]
if mdpOptions.get('implicit_solvent') == 'gbsa':
    system = top.createSystem(nonbondedMethod=NoCutoff, implicitSolvent=OBC2)
elif mdpOptions.get('coulombtype') == 'pme':
    system = top.createSystem(nonbondedMethod=PME, ewaldErrorTolerance=5e-5, nonbondedCutoff=0.9*nanometers)
else:
    system = top.createSystem(nonbondedMethod=NoCutoff)
for i,f in enumerate(system.getForces()):
    f.setForceGroup(i)
integrator = VerletIntegrator(0.001*picoseconds)
simulation = Simulation(top.topology, system, integrator, Platform.getPlatformByName(platformName), properties)
simulation.context.setPositions(gro.positions)
simulation.context.applyConstraints(1e-6)
state = simulation.context.getState(getForces=True)
forces = state.getForces().value_in_unit(kilojoules_per_mole/nanometer)

# Report results.

relativeDiff = [2*norm(f1-f2)/(norm(f1)+norm(f2)) for f1, f2 in zip(forces, gromacsForces)]
relativeDiff.sort()
print('Relative force error')
print('max:', relativeDiff[-1])
print('mean: ', sum(relativeDiff)/len(relativeDiff))
print('median:', relativeDiff[len(relativeDiff)/2])
projection = [dot(f1, f2)/dot(f1, f1) for f1, f2 in zip(forces, gromacsForces)]
print()
print('Projection')
print('mean:', sum(projection)/len(projection))
print('min:', min(projection))
print('max:', max(projection))
print()
print('OpenMM component energies')
for i,f in enumerate(system.getForces()):
    state = simulation.context.getState(getEnergy=True, groups={i})
    print(f.__class__, state.getPotentialEnergy())
print()
print('Gromacs component energies')
for term in termEnergies:
    print(term, termEnergies[term])

