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
import os
from pathlib import Path

import ase
from ase import io
from ase.visualize import view

import matplotlib.pyplot as plt

from tqdm import tqdm

def dihedral3(p):
    b = p[:-1] - p[1:]
    b[0] *= -1
    v = np.array( [np.cross(v,b[1]) for v in [b[0], b[2]] ] )
    # Normalize vectors
    v /= np.sqrt(np.einsum('...i,...i', v, v)).reshape(-1,1)
    return np.degrees(np.arccos( v[0].dot(v[1]) ))

def cond(x):
    return 6 if x > 10 else  x

def vis_pdb(pdb):
    atomic_numbers = pdb.get_atomic_numbers()
    print(atomic_numbers)
    #  hack to get around the color scheme and atom type names
    #. eg. CYA for alpha carbon types
    pdb.set_atomic_numbers([ cond(_)
                            for _ in atomic_numbers])
    return view(pdb, viewer="x3d")

def mdcalc(namepdb, s, t, steps, timesteps, n, i):
    wb=False
    print(f"{s+str(t)+str(wb)+str(n)}")
    
    
    initial_pdb = io.read(f'/home/himmelreich/PCProject/charmmdynamics/analysispdb/{s+str(t)+str(wb)+str(n)}/{namepdb}')

    OUTPUT = "opt_min.xyz"
    atoms = initial_pdb
    calc = NNCalculator(
        checkpoint="models/Ac-Ala3-NHMe-coff10-elec-disp-wf53-cw14-alldata-b", #load the model you want to used
        atoms=atoms,
        charge=0,
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

    #optimize
    opt = BFGS(atoms)
    #opt.run(10**6)

    #write output file
    #write(OUTPUT, atoms)


    TEMP = t[0]
    CHARGE = 0
    STEPS = steps
    TIMESTEP = timesteps
    FRICTION = 0.02
    INTERVAL = 20
    filename = f"{s}_dynamics_{i}"

    #run an optimization
    #BFGS(atoms).run(fmax=0.01)

    # Set the momenta corresponding to a temperature T
    MaxwellBoltzmannDistribution(atoms, TEMP * units.kB)
    ZeroRotation(atoms)
    Stationary(atoms)

    # define the algorithm for MD: here Langevin alg. with with a time step of 0.1 fs,
    # the temperature T and the friction coefficient to 0.02 atomic units.
    dyn = Langevin(atoms, TIMESTEP * units.fs, TEMP * units.kB, FRICTION)
    #dyn = VelocityVerlet(atoms, args.timestep * units.fs)


    def printenergy(a=atoms):  # store a reference to atoms in the definition.
        """Function to print the potential, kinetic and total energy."""
        epot = a.get_potential_energy() / len(a)
        ekin = a.get_kinetic_energy() / len(a)
        print('Energy per atom: Epot = %.3feV  Ekin = %.3feV (T=%3.0fK)  '
              'Etot = %.3feV' % (epot, ekin, ekin / (1.5 * units.kB), epot + ekin))
    
    
    if os.path.exists("/home/himmelreich/PCProject/Ac-Ala3-NHMe-MD/trajs"):
        pass
    else:
        os.mkdir("/home/himmelreich/PCProject/Ac-Ala3-NHMe-MD/trajs")

    # save the positions of all atoms after every Xth time step.
    traj = Trajectory(f"trajs/{str(TEMP)+ 'K_md_' + filename + '.traj'}", 'w', atoms)

    """
    #equilibration
    for i in range(10000):
        if i%100 == 0:
            print("Equilibration Step: ", i)
        dyn.run(1)

    """
    start_time = time.time()
    # run the dynamics
    for i in tqdm(range(STEPS)):
        dyn.run(1)
        if i%INTERVAL == 0:
            #epot = atoms.get_potential_energy() / len(atoms)
            #ekin = atoms.get_kinetic_energy() / len(atoms)
            # print("Production Step: ", i)
            traj.write()

    end_time = time.time()


    # Calculate elapsed time
    elapsed_time = end_time - start_time
    print("Elapsed time: ", elapsed_time)
    print("Time per step: ", elapsed_time / STEPS)

    #intraj = Trajectory("298K_md_dynamics.traj")

def get_list_of_pdb(s, t, n):
    wb=False
    charmdynhome = Path("/home/himmelreich/PCProject/charmmdynamics")
    pdb = charmdynhome / f"analysispdb/{s+str(t)+str(wb)+str(n)}"
    list_of_pdbs = os.listdir(pdb)
    return list_of_pdbs
    
    
    