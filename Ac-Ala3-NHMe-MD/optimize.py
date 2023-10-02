#!/usr/bin/env python3
import argparse
from ase import Atoms
from ase.io import read, write
import argparse
import numpy as np
from ase.optimize import *
from NNCalculator.NNCalculator import *
from ase.optimize import *

#parse command line arguments
parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser._action_groups.pop()
required = parser.add_argument_group("required arguments")
required.add_argument("-i", "--input",   type=str,   help="input xyz",  required=True)
required.add_argument("-o", "--output",  type=str,   help="output xyz", required=True)
optional = parser.add_argument_group("optional arguments")
optional.add_argument("--charge",  type=float, help="total charge",            default=0.0)
optional.add_argument("--fmax",    type=float, help="optimizer tolerance",     default=0.04)
args = parser.parse_args()
print("input ", args.input)
print("output", args.output)

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

#optimize
opt = BFGS(atoms)
opt.run(fmax=args.fmax)

#write output file
write(args.output, atoms)
