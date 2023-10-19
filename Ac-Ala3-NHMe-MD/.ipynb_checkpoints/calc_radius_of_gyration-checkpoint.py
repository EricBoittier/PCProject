#!/usr/bin/env python3

# imports
import argparse
from ase import Atoms
from ase.visualize import view
from ase.io import read, write
from ase.calculators.mopac import MOPAC
from ase.optimize import *
from ase import units
from ase.md.velocitydistribution import MaxwellBoltzmannDistribution
from ase.md.langevin import Langevin
from ase.io.trajectory import Trajectory
from os.path import splitext
import io, os
import numpy as np

#parse command line arguments
parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser._action_groups.pop()
required = parser.add_argument_group("required arguments")
required.add_argument("-i", "--input",   type=str,   help="input xyz",  required=True)

    
def radgyr(atoms, masses, total_mass=None):
    # coordinates change for each frame
    coordinates = atoms.get_positions()
    center_of_mass = atoms.get_center_of_mass()

    # get squared distance from center
    ri_sq = (coordinates-center_of_mass)**2
    # sum the unweighted positions
    sq = np.sum(ri_sq, axis=1)
    sq_x = np.sum(ri_sq[:,[1,2]], axis=1) # sum over y and z
    sq_y = np.sum(ri_sq[:,[0,2]], axis=1) # sum over x and z
    sq_z = np.sum(ri_sq[:,[0,1]], axis=1) # sum over x and y

    # make into array
    sq_rs = np.array([sq, sq_x, sq_y, sq_z])

    # weight positions
    rog_sq = np.sum(masses*sq_rs, axis=1)/total_mass

    # square root and return
    return np.sqrt(rog_sq)





args = parser.parse_args()
filename, extension = splitext(args.input)

traj = Trajectory(args.input)



atoms = traj[0]

del atoms[[atom.symbol == 'H' for atom in atoms]]
masses = atoms.get_masses()
mtot = np.sum(masses)

rg = np.zeros((len(traj), 4))
for ii in range(len(traj)):
    atoms = traj[ii]

    del atoms[[atom.symbol == 'H' for atom in atoms]]

    rg[ii, :] = radgyr(atoms, masses, mtot)
    
import matplotlib.pyplot as plt
fig,ax = plt.subplots(figsize=(8,8))

# Fontsize
SMALL_SIZE = 14
MEDIUM_SIZE = 16
BIGGER_SIZE = 18

#plt.rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
plt.rc('font', size=SMALL_SIZE, weight='bold')  # controls default text sizes
plt.rc('axes', titlesize=SMALL_SIZE)     # fontsize of the axes title
plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('legend', fontsize=SMALL_SIZE)    # legend fontsize
plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title

plt.ylabel(r'$\mathrm{Radius\;\ of\;\ Gyration} \;\, \mathrm{(\AA)}$', fontsize=12)
plt.xlabel(r'$\mathrm{Time} \;\, \mathrm{(ns)}$', fontsize=12)

#colors = ["#003f5c", "#2f4b7c", "#665191", "#a05195", "#d45087", "#f95d6a", "#ff7c43", "#ffa600"]
colors = ["#003f5c", "#665191", "#d45087", "#ff7c43"]

xaxis = np.linspace(0, 1, 100)
plt.plot(xaxis, rg[:, 0], color=colors[0])

plt.xlim(0.0, 0.7)
plt.ylim(2.0, 4.0)
plt.savefig("Radius_gyration.png",bbox_inches='tight',dpi=300)

print(rg)

plt.show()

