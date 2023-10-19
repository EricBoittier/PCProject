import pandas as pd
import os
import numpy as np
import matplotlib.pyplot as plt
import ase
from ase import io
from ase.visualize import view
from pathlib import Path
import matplotlib.pyplot as plt

from prody import *

from sklearn.neighbors import KernelDensity


''' functions we need:
    
'''

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

def sillyfunction():
    return "silly"

def dihedral3(p):
    b = p[:-1] - p[1:]
    b[0] *= -1
    v = np.array( [np.cross(v,b[1]) for v in [b[0], b[2]] ] )
    # Normalize vectors
    v /= np.sqrt(np.einsum('...i,...i', v, v)).reshape(-1,1)
    return np.degrees(np.arccos( v[0].dot(v[1]) ))


def new_dihedral(p):
    """Praxeolitic formula
    1 sqrt, 1 cross product"""
    p0 = p[0]
    p1 = p[1]
    p2 = p[2]
    p3 = p[3]

    b0 = -1.0*(p1 - p0)
    b1 = p2 - p1
    b2 = p3 - p2

    # normalize b1 so that it does not influence magnitude of vector
    # rejections that come next
    b1 /= np.linalg.norm(b1)

    # vector rejections
    # v = projection of b0 onto plane perpendicular to b1
    #   = b0 minus component that aligns with b1
    # w = projection of b2 onto plane perpendicular to b1
    #   = b2 minus component that aligns with b1
    v = b0 - np.dot(b0, b1)*b1
    w = b2 - np.dot(b2, b1)*b1

    # angle between v and w in a plane is the torsion angle
    # v and w may not be normalized but that's fine since tan is y/x
    x = np.dot(v, w)
    y = np.dot(np.cross(b1, v), w)
    return np.degrees(np.arctan2(y, x))


dihedral3 = new_dihedral

# path the data
charmdynhome = Path("/home/himmelreich/PCProject/charmmdynamics")
#  dictionary to link the aminoacid to the atom indices used to calculate phi and psi
seq_to_atomID = {"ala": ([4,6,8,14], [6,8,14,16], [8, 14, 16, 18]), 
                 "gly": ([4,6,8,11], [6,8,11,13], [8, 11, 13, 15]), 
                 "ser": ([4,6,8,15], [6,8,15,17], [8, 15, 17, 19]), 
                 "thr": ([4,6,8,18], [6,8,18,20], [8, 18, 20, 22])}


def get_stats_from_key(key, s):
    """ s = 'ala' or 'ser' 
    """
    # name of the dcdfile (/ joins the path)
    #  f"{x}" format string
    dcd = charmdynhome / f"dcd/name_{key}.dcd"
    # name of the pdbfile
    pdb = charmdynhome / f"pdb/{s.lower()}-{key}-end.pdb"
    
    phi_indx, psi_indx, omega_index = seq_to_atomID[s]
    
    dcd = str(dcd)
    pdb = str(pdb)
    structure = parsePDB(pdb)
    traj = Trajectory(dcd)
    traj.link(structure)
    
    nter = structure.select('name NT and resnum 1')
    cter = structure.select('name CAY and resnum 1')
    protein = structure.select('resnum 1')
    
    rgyr = np.zeros(traj.numFrames())
    e2e = np.zeros(traj.numFrames())
    psi = np.zeros(traj.numFrames())
    coords = []
    
    for i, frame in enumerate(traj):
        e2e[i] = calcDistance(nter, cter)
        rgyr[i] = calcGyradius(protein)
        coords.append(getCoords(protein))
        
    phidihs = [dihedral3(coords[i][phi_indx]) for i
        in range(len(coords))]
    psidihs = [dihedral3(coords[i][psi_indx]) for i
        in range(len(coords))]
    omegadihs = [dihedral3(coords[i][omega_index]) for i
        in range(len(coords))]
    
    ramachandran = pd.DataFrame({
        "sequence":s,
        "key": key,
        "phi":phidihs, 
        "psi": psidihs,
        "omega": omegadihs,               
        "e2e": e2e, 
        "rgyr": rgyr})
    ramachandran.to_csv(f"datacsv/{key}.csv", index=False)
    return ramachandran
    
    
             
def load_csv(s, t, wb):
    csvpath = charmdynhome / "csv"
    csvfiles = csvpath.glob("*.csv")
    dcdpath = charmdynhome / "dcd"
    dcdfiles = dcdpath.glob("*.dcd")
    resultsdf = pd.concat([pd.read_csv(_) for _ in list(csvfiles)])
    resultsdf["time (ps)"] = resultsdf["dt"] * resultsdf["n"]
    selectedcsv = resultsdf[(resultsdf["s"] == s.upper()) & (resultsdf["wb"] == wb) & (resultsdf["t"].isin(t))]
    print("works")
    return selectedcsv

def printgraph(s,t,wb):
    dataframe = load_csv(s, t, wb)
    dataframes = []
    #print(dataframe)
    for row in dataframe.iterrows():

        key = row[1][0]
        s = row[1][1].lower()
        tmpdf = get_stats_from_key(key,s)
        dataframes.append(tmpdf)
    master = pd.concat(dataframes)
    #print(master)
    
    if wb:
        wbtitle="_wb"
    else:
        wbtitle=""
    plt.scatter(master["phi"], master["psi"], alpha=0.01, color='tab:blue')
    plt.title(f"{s.upper()}_{str(t)}K{wbtitle}")
    plt.xlabel("$\phi$ Angle ($^{\circ}$)")
    plt.ylabel("$\psi$ Angle ($^{\circ}$)")
    plt.xlim(-180,180)
    plt.ylim(-180,180)
    plt.savefig(f"plot/{s}_{str(t)}K_{wbtitle}_title.png", dpi=300)
    return master

def getstartingpdb (s, t, wb, m, inverse=False):
    datacsvpath = charmdynhome / "datacsv"
    allramachandran = datacsvpath.glob("*.csv")
    keys = load_csv(s, t, wb)["Unnamed: 0"]
    
    listoffiles = []
    for key in keys:
        #print(key)
        csvfile=datacsvpath/f"{key}.csv"
        listoffiles.append(csvfile)
    # combine the csvs
    selected = pd.concat([pd.read_csv(_) for _ in listoffiles])
    
    # sample m samples randomly
    if not inverse:
        frames =selected.sample(n = m)
        addtitle = ""
    else:
        # sample inverse prob.
        X = selected[["phi", "psi"]]
        kde = KernelDensity(kernel='gaussian', bandwidth=0.2).fit(X)
        selected["invloglik"] = abs( 1 / kde.score_samples(X))
        #frames = selected.sample(n = m, weights=invprobs)
        frames = selected[probabilities['invloglik'] < 0.2].sample(n=m)
        addtitle = "inverse-"
        
        
    
    
    
    
    new_dir_path = f"analysispdb/{s+str(t)+str(wb)+str(m)}"
    if os.path.exists(new_dir_path):
        pass
    else:
        os.mkdir(new_dir_path)
    for _ in frames.iterrows():
        #print(_,"lol")
        key=_[1][1]
        print(key)
        frame=_[0]
        #print(frame)
        
        # name of the dcdfile (/ joins the path)
        #  f"{x}" format string
        dcd = charmdynhome / f"dcd/name_{key}.dcd"
        # name of the pdbfile
        pdb = charmdynhome / f"pdb/{s.lower()}-{key}-end.pdb"


        dcd = str(dcd)
        pdb = str(pdb)
        structure = parsePDB(pdb)
        traj = Trajectory(dcd)
        traj.link(structure)
        writePDB(f"analysispdb/{s+str(t)+str(wb)+str(m)}/{key}-{addtitle}frame{frame}.pdb", traj[frame].getAtoms())
    return frames
        
    
    








