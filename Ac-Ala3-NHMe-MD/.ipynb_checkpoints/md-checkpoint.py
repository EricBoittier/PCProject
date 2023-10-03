#!/usr/bin/env python3

# imports
import argparse
from ase import Atoms
from ase.io import read, write
from ase.optimize import *
from ase import units
from ase.md.velocitydistribution import MaxwellBoltzmannDistribution, ZeroRotation, Stationary
from ase.md.langevin import Langevin
from ase.io.trajectory import Trajectory
from os.path import splitext
from ase.visualize import view
from ase.md.verlet import VelocityVerlet
import numpy as np
from NNCalculator.NNCalculator import *
import time


#parse command line arguments
parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser._action_groups.pop()
required = parser.add_argument_group("required arguments")
required.add_argument("-i", "--input",   type=str,   help="input xyz",  required=True)
optional = parser.add_argument_group("optional arguments")
optional.add_argument("--charge",  type=float, help="total charge", default=0.0)
optional.add_argument("--temperature", type=float, help="Set momenta corresponding to a temperature T", default=300)
optional.add_argument("--timestep",  type=float, help="timestep for Langevin algorithm", default=0.5)
optional.add_argument("--friction",  type=float, help="friction coeff for Langevin algorithm", default=1e-3)
optional.add_argument("--steps",  type=int, help="number of steps in md", default=2000000)
optional.add_argument("--interval",  type=float, help="interval for saving snapshots", default=200)

args = parser.parse_args()
filename, extension = splitext(args.input)


#read input file 
atoms = read(args.input)

#setup calculator object, which in this case is the NN calculator
#it is important that it is setup with the settings as used in the
#training procedure.
#setup calculator
calc = NNCalculator(
    checkpoint="models/Ac-Ala3-NHMe-coff10-elec-disp-wf53-cw14-alldata-b", #load the model you want to used
    atoms=atoms,
    charge=args.charge,
    F=128,
    K=64,
    num_blocks=5,
    num_residual_atomic=2,
    num_residual_interaction=3,
    num_residual_output=1,
    sr_cut=10.0,
    use_electrostatic=True,
    use_dispersion=True,
    s6=1.0000,                    #s6 coefficient for d3 dispersion, by default is learned
    s8=2.3550,                    #s8 coefficient for d3 dispersion, by default is learned
    a1=0.5238,                    #a1 coefficient for d3 dispersion, by default is learned
    a2=3.5016)                   #a2 coefficient for d3 dispersion, by default is learned)

#attach the calculator object to the atoms object
atoms.set_calculator(calc)               

#run an optimization
BFGS(atoms).run(fmax=0.01)


# Set the momenta corresponding to a temperature T
MaxwellBoltzmannDistribution(atoms, args.temperature * units.kB)
ZeroRotation(atoms)
Stationary(atoms)

# define the algorithm for MD: here Langevin alg. with with a time step of 0.1 fs,
# the temperature T and the friction coefficient to 0.02 atomic units.
dyn = Langevin(atoms, args.timestep * units.fs, args.temperature * units.kB, args.friction)
#dyn = VelocityVerlet(atoms, args.timestep * units.fs)


def printenergy(a=atoms):  # store a reference to atoms in the definition.
    """Function to print the potential, kinetic and total energy."""
    epot = a.get_potential_energy() / len(a)
    ekin = a.get_kinetic_energy() / len(a)
    print('Energy per atom: Epot = %.3feV  Ekin = %.3feV (T=%3.0fK)  '
          'Etot = %.3feV' % (epot, ekin, ekin / (1.5 * units.kB), epot + ekin))
    
# save the positions of all atoms after every Xth time step.
traj = Trajectory(str(args.temperature)+ 'K_md_' + filename + '.traj', 'w', atoms)

"""
#equilibration
for i in range(10000):
    if i%100 == 0:
        print("Equilibration Step: ", i)
    dyn.run(1)

"""
start_time = time.time()
# run the dynamics
for i in range(args.steps):
    dyn.run(1)
    if i%args.interval == 0:
        #epot = atoms.get_potential_energy() / len(atoms)
        #ekin = atoms.get_kinetic_energy() / len(atoms)
        print("Production Step: ", i)
        traj.write()

end_time = time.time()


# Calculate elapsed time
elapsed_time = end_time - start_time
print("Elapsed time: ", elapsed_time)
print("Time per step: ", elapsed_time / args.steps)


