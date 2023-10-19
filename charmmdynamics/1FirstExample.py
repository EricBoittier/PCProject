import os
import sys
import pandas as pd

# These are a subset of the pycharmm modules that were installed when
# pycharmm was installed in your python environment
import pycharmm
import pycharmm.generate as gen
import pycharmm.ic as ic
import pycharmm.coor as coor
import pycharmm.energy as energy
import pycharmm.dynamics as dyn
import pycharmm.nbonds as nbonds
import pycharmm.minimize as minimize
import pycharmm.crystal as crystal
import pycharmm.image as image
import pycharmm.psf as psf
import pycharmm.read as read
import pycharmm.write as write
import pycharmm.settings as settings
import pycharmm.cons_harm as cons_harm
import pycharmm.cons_fix as cons_fix
import pycharmm.select as select
import pycharmm.shake as shake
from pycharmm.lib import charmm as libcharmm


import argparse

import uuid

from dyna_dicts import get_dynamics_dict


#  settings
cutnb = 12
cutim = cutnb
ctofnb = cutnb - 1.0
ctonnb = cutnb - 3.0

pycharmm.lingo.charmm_script("bomlev -9")


def run_md(
    pdbid="name",
    useomm=False,
    useblade=False,
    nequil=10000,
    nsteps=5000,
    nsavc=100,
    leap=True,
    lang=False,
    DT=0.002,
    T=298.0,
    uuid=0,
):
    if useomm:
        append = "omm"
    elif useblade:
        append = "blade"

    # dyn.set_fbetas(np.full((psf.get_natom()),1.0,dtype=float))

    res_file = pycharmm.CharmmFile(
        file_name=f"res/{pdbid}.res", file_unit=2, formatted=True, read_only=False
    )
    dcd_file = pycharmm.CharmmFile(
        file_name=f"res/{pdbid}.dcd".format(pdbid),
        file_unit=3,
        formatted=False,
        read_only=False,
    )

    dyndict = get_dynamics_dict(res_file=res_file, dcd_file=dcd_file, dynatype="heat", timestep=DT, temp=T, nstep=nequil )

    my_dyn = pycharmm.DynamicsScript(**dyndict)

    my_dyn.run()

    res_file.close()
    dcd_file.close()
    # open unit 2 write form name res/{}.res
    res_file = pycharmm.CharmmFile(
        file_name=f"res/{pdbid}.res", file_unit=2, formatted=True, read_only=False
    )
    lam_file = pycharmm.CharmmFile(
        file_name=f"res/{pdbid}.lam", file_unit=3, formatted=False, read_only=False
    )
    # open unit 1 write file name dcd/{}.dcd
    dcd_file = pycharmm.CharmmFile(
        file_name=f"dcd/{pdbid}_{uuid}.dcd",
        file_unit=1,
        formatted=False,
        read_only=False,
    )

    my_dyn = pycharmm.DynamicsScript(
        leap=leap,
        lang=lang,
        start=False,
        restart=True,
        nstep=nsteps,
        timest=DT,
        firstt=T,
        finalt=T,
        tbath=T,
        tstruc=T,
        teminc=0.0,
        twindh=0.0,
        twindl=0.0,
        iunwri=res_file.file_unit,
        iunrea=res_file.file_unit,
        iuncrd=dcd_file.file_unit,
        # iunlam=lam_file.file_unit,
        inbfrq=-1,
        imgfrq=-1,
        iasors=0,
        iasvel=1,
        ichecw=0,
        iscale=0,
        iscvel=0,
        echeck=-1,
        nsavc=nsavc,
        nsavv=0,
        nsavl=0,
        ntrfrq=0,
        isvfrq=nsavc,
        iprfrq=2 * nsavc,
        nprint=nsavc,
        ihtfrq=0,
        ieqfrq=0,
        ilbfrq=0,
        ihbfrq=0,
        omm=useomm,
        blade=useblade,
    )
    my_dyn.run()

    res_file.close()
    lam_file.close()
    dcd_file.close()


if __name__ == "__main__":

    uid = uuid.uuid4()

    par = argparse.ArgumentParser()
    par.add_argument("--s", "-s", default="ALA")
    par.add_argument("--t", "-t", default=298.0)
    par.add_argument("--n", "-n", default=1000000)
    par.add_argument("--o", "-o", default="output")
    par.add_argument("--dt", "-dt", default=0.002)
    par.add_argument("--wb", "--wb", default=False)
    args = par.parse_args()

    pd.DataFrame(vars(args), index=[uid]).to_csv(f"csv/{uid}.csv")

    # Read in the topology (rtf) and parameter file (prm) for proteins
    # equivalent to the CHARMM scripting command: read rtf card name toppar/top_all36_prot.rtf
    read.rtf("../toppar/top_all36_prot.rtf")
    # equivalent to the CHARMM scripting command: read param card flexible name toppar/par_all36m_prot.prm
    read.prm("../toppar/par_all36m_prot.prm", flex=True)

    # stream in the water/ions parameter using the pycharmm.lingo module
    # equivalent to the CHARMM scripting command: stream toppar/toppar_water_ions.str
    pycharmm.lingo.charmm_script("stream ../toppar/toppar_water_ions.str")
    # end toppar/toppar_water_ions.str

    read.sequence_string(args.s)

    # equivalent to the CHARMM scripting command: generate ADP first ACE last CT3 setup
    gen.new_segment(seg_name="ADP", first_patch="ACE", last_patch="CT3", setup_ic=True)

    # equivalent to the CHARMM scripting command: ic param
    ic.prm_fill(replace_all=False)
    # equivalent to the CHARMM scripting command: ic seed 1 CAY 1 CY 1 N
    ic.seed(res1=1, atom1="CAY", res2=1, atom2="CY", res3=1, atom3="N")
    # equivalent to the CHARMM scripting command: ic build
    ic.build()

    # The coor orie command is useful to expose since it allows one to
    # orient the system in preparation for other calculations
    # equivalent to the CHARMM scripting command: coor orient
    coor.orient(by_rms=False, by_mass=False, by_noro=False)
    # equivalent to the CHARMM scripting command: print coor
    coor.show()
    # If pdb directory doesn't alrady exist make it here.
    if not os.path.isdir("pdb"):
        os.system("mkdir pdb")
    # equivalent to the CHARMM scripting command: write coor pdb name pdb/initial.pdb

    write.coor_pdb(f"pdb/initial-{args.o}-{uid}.pdb")

    my_nbonds = pycharmm.NonBondedScript(
        cutnb=18.0,
        ctonnb=13.0,
        ctofnb=17.0,
        eps=1.0,
        cdie=True,
        atom=True,
        vatom=True,
        fswitch=True,
        vfswitch=True,
    )

    # Implement these non-bonded parameters by "running" them.
    my_nbonds.run()

    # equivalent CHARMM scripting command: minimize abnr nstep 1000 tole 1e-3 tolgr 1e-3
    minimize.run_abnr(nstep=1000, tolenr=1e-3, tolgrd=1e-3)
    # equivalent CHARMM scripting command: energy
    energy.show()

    # write coor pdb name pdb/adp.pdb
    write.coor_pdb(f"pdb/{args.s.lower()}-{uid}.pdb")
    # write psf card name pdb/adp.psf
    write.psf_card(f"pdb/{args.s.lower()}-{uid}.psf")

    # delete atom select all end
    psf.delete_atoms(pycharmm.SelectAtoms().all_atoms())

    # read psf card name pdb/adp.psf
    read.psf_card(f"pdb/{args.s.lower()}-{uid}.psf")

    # read coor pdb name pdb/adp.pdb resid
    adp_pdb_file = f"pdb/{args.s.lower()}-{uid}.pdb"
    read.pdb(adp_pdb_file, resid=True)

    # read psf card name pdb/adp.psf
    read.psf_card(f"pdb/{args.s.lower()}-{uid}.psf")
    adp_pdb_file = f"pdb/{args.s.lower()}-{uid}.pdb"
    read.pdb(adp_pdb_file, resid=True)

    stats = coor.stat()
    print(stats)

    if args.wb:
        CONVPDB = "/home/himmelreich/PCProject/toolset/perl/convpdb.pl"

        # CHARMM scripting command: system "convpdb.pl -solvate -cutoff 10 -cubic -out charmm22 pdb/adp.pdb
        # | convpdb.pl -segnames adp ala -nsel TIP3 > pdb/wt00.pdb"
        solvate_command = (
            f"{CONVPDB} -solvate -cutoff 10 -cubic -out charmm22 {adp_pdb_file} | "
        )
        solvate_command += f"{CONVPDB} -segnames -nsel TIP3 > pdb/wt00.pdb"
        # run the command as a system subprocess
        os.system(solvate_command)

        #  how many water molecules?
        N = int(open("pdb/wt00.pdb").readlines()[-3].split()[4]) - 1 # minus the peptide

        read.sequence_string("TIP3 " * N)
        # read.sequence_pdb('water_box.pdb')
        gen.new_segment("WT00", angle=False, dihedral=False)
        read.pdb("pdb/wt00.pdb")

        stats = coor.stat()
        print(stats)

        # boxsize
        xsize = stats["xmax"] - stats["xmin"]
        ysize = stats["ymax"] - stats["ymin"]
        zsize = stats["zmax"] - stats["zmin"]
        boxsize = max(xsize, ysize, zsize)

        # half box size
        boxhalf = boxsize / 2.0

        # CHARMM scripting: crystal define cubic @boxsize @boxsize @boxsize 90 90 90
        crystal.define_cubic(boxsize)
        # CHARMM scripting: crystal build cutoff @boxhalf noper 0
        crystal.build(boxhalf)

        # Turn on image centering - bysegment for peptide, by residue for solvent
        # CHARMM scripting: image byseg xcen 0 ycen 0 zcen 0 select segid adp end
        image.setup_segment(0.0, 0.0, 0.0, "ADP")
        # CHARMM scripting: image byres xcen 0 ycen 0 zcen 0 select resname tip3 end
        image.setup_residue(0.0, 0.0, 0.0, "TIP3")

        # Now specify nonbonded cutoffs for solvated box
        cutnb = min(boxhalf, 12)
        cutim = cutnb
        ctofnb = cutnb - 1.0
        ctonnb = cutnb - 3.0

    # Another nbonds example
    # CHARMM scripting: nbonds cutnb @cutnb cutim @cutim ctofnb @ctofnb ctonnb @ctonnb -
    #        inbfrq -1 imgfrq -1
    pycharmm.NonBondedScript(
        cutnb=cutnb,
        cutim=cutim,
        ctonnb=ctonnb,
        ctofnb=ctofnb,
        eps=1.0,
        cdie=True,
        atom=True,
        vatom=True,
        fswitch=True,
        vfswitch=True,
        inbfrq=-1,
        imgfrq=-1,
    ).run()

    # Fix the peptide and minimize the solvent to "fit"
    # CHARMM scripting: cons fix select segid adp end
    cons_fix.setup(pycharmm.SelectAtoms(seg_id="ADP"))

    # Minimize the solvent positions with periodic boundary conditions using steepest descents
    # CHARMM scripting: mini sd nstep 200 tole 1e-3 tolgrd 1e-3
    minimize.run_sd(nstep=200, tolenr=1e-3, tolgrd=1e-3)

    # Turn off fixed atoms
    # CHARMM scripting: cons fix select none end
    cons_fix.turn_off()

    # Write the psf and coordinates for the solvated peptide
    # write psf card name pdb/adp+wat.psf
    write.psf_card(f"pdb/{args.s.lower()}-{uid}-min.psf")
    # write coor pdb name pdb/adp+wat_min.pdb
    write.coor_pdb(f"pdb/{args.s.lower()}-{uid}-min.pdb")

    run_md(uuid=uid, nsteps=args.n, T=args.t, DT=args.dt)
    
    # Write the psf and coordinates for the solvated peptide
    # write psf card name pdb/adp+wat.psf
    write.psf_card(f"pdb/{args.s.lower()}-{uid}-end.psf")
    # write coor pdb name pdb/adp+wat_min.pdb
    write.coor_pdb(f"pdb/{args.s.lower()}-{uid}-end.pdb")
