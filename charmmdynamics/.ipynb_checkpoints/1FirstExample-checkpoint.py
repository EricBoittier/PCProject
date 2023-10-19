import os
import sys

# This should specify the path to the install pyCHARMM library
# pyCHARMM_LIB = '/Users/brookscl/charmm/c47-dev-release/install-pycharmm-nompi'
# if os.getenv('CHARMM_LIB_DIR') == None:
#     os.environ['CHARMM_LIB_DIR'] = pyCHARMM_LIB
#     print(os.getenv('CHARMM_LIB_DIR'))

os.environ["PYCHARMMPATH"] = "/home/boittier/dev-release-c48a-fmdcm-hotfix"
os.environ["PYTHONPATH"] = "/home/boittier/dev-release-c48a-fmdcm-hotfix/tool/pycharmm:" 
os.environ["CHARMM_LIB_DIR"] = "/home/boittier/dev-release-c48a-fmdcm-hotfix/lib/"
os.environ["CHARMM_DATA_DIR"] = "/home/boittier/dev-release-c48a-fmdcm-hotfix/test/data/"

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

# Read in the topology (rtf) and parameter file (prm) for proteins
# equivalent to the CHARMM scripting command: read rtf card name toppar/top_all36_prot.rtf
read.rtf('../toppar/top_all36_prot.rtf')
# equivalent to the CHARMM scripting command: read param card flexible name toppar/par_all36m_prot.prm
read.prm('../toppar/par_all36m_prot.prm', flex=True)

# stream in the water/ions parameter using the pycharmm.lingo module
# equivalent to the CHARMM scripting command: stream toppar/toppar_water_ions.str
pycharmm.lingo.charmm_script('stream ../toppar/toppar_water_ions.str')
# end toppar/toppar_water_ions.str


read.sequence_string('ALA')

# equivalent to the CHARMM scripting command: generate ADP first ACE last CT3 setup
gen.new_segment(seg_name='ADP', first_patch='ACE', last_patch='CT3', setup_ic=True)

# equivalent to the CHARMM scripting command: ic param
ic.prm_fill(replace_all=False)
# equivalent to the CHARMM scripting command: ic seed 1 CAY 1 CY 1 N
ic.seed(res1=1, atom1='CAY', res2=1, atom2='CY', res3=1, atom3='N')
# equivalent to the CHARMM scripting command: ic build
ic.build()

# The coor orie command is useful to expose since it allows one to
# orient the system in preparation for other calculations
# equivalent to the CHARMM scripting command: coor orient
coor.orient(by_rms=False,by_mass=False,by_noro=False)
# equivalent to the CHARMM scripting command: print coor
coor.show()
# If pdb directory doesn't alrady exist make it here.
if not os.path.isdir('pdb'): os.system('mkdir pdb')
# equivalent to the CHARMM scripting command: write coor pdb name pdb/initial.pdb
write.coor_pdb('pdb/initial.pdb')

my_nbonds = pycharmm.NonBondedScript(
    cutnb=18.0, ctonnb=13.0, ctofnb=17.0,
    eps=1.0,
    cdie=True,
    atom=True, vatom=True,
    fswitch=True, vfswitch=True)

# Implement these non-bonded parameters by "running" them.
my_nbonds.run()

# equivalent CHARMM scripting command: minimize abnr nstep 1000 tole 1e-3 tolgr 1e-3
minimize.run_abnr(nstep=1000, tolenr=1e-3, tolgrd=1e-3)
# equivalent CHARMM scripting command: energy
energy.show()


# write coor pdb name pdb/adp.pdb
write.coor_pdb('pdb/adp.pdb')
# write psf card name pdb/adp.psf
write.psf_card('pdb/adp.psf')

# delete atom select all end
psf.delete_atoms(pycharmm.SelectAtoms().all_atoms())

# read psf card name pdb/adp.psf
read.psf_card('pdb/adp.psf')

# read coor pdb name pdb/adp.pdb resid
adp_pdb_file = 'pdb/adp.pdb'
read.pdb(adp_pdb_file, resid=True)

print(coor.stat())


PACKMOL = "/home/boittier/PCProject/packmol/packmol"
PDBPATH = "pdb/adp.pdb"

job_file = "water_box.inp"
job_str = f"""#
# A box with water
#

tolerance 2.0
filetype pdb
output water_box.pdb

#structure {PDBPATH}
#  center
#  fixed 10. 10. 10. 0. 0. 0.
#end structure

structure water.pdb
  number 1000
  inside box 0. 0. 0. 40. 40. 40.
end structure

"""
with open(job_file, "w") as f:
    f.write(job_str)
    
    
os.system(f"{PACKMOL} < water_box.inp;" +
          "sed -i 's/HOH /TIP3/g' water_box.pdb;"+
         "sed -i 's/HETATM/ATOM/g' water_box.pdb"
         )

stats = coor.stat()
print(stats)  

psf.delete_atoms(pycharmm.SelectAtoms().all_atoms())




# stats = coor.stat()
# print(stats)

# read sequ pdb name pdb/wt00.pdb
# read.sequence_pdb('water_box.pdb', resid="TIP3")
# read psf card name pdb/adp.psf
read.psf_card('pdb/adp.psf')
adp_pdb_file = 'pdb/adp.pdb'
read.pdb(adp_pdb_file, resid=True)

stats = coor.stat()
print(stats)

read.sequence_string( "TIP3 " * 1000)
# read.sequence_pdb('water_box.pdb')
gen.new_segment('WT00', angle=False, dihedral=False)
read.pdb('water_box.pdb')

#Another example of the generate command
#generate wt00 noangle nodihedral

# read psf card name pdb/adp.psf
# read.psf_card('pdb/adp.psf')

# read coor pdb name pdb/adp.pdb resid
# adp_pdb_file = 'pdb/adp.pdb'
# read.pdb(adp_pdb_file, resid=True)


stats = coor.stat()
print(stats)

# boxsize
xsize = stats['xmax'] - stats['xmin']
ysize = stats['ymax'] - stats['ymin']
zsize = stats['zmax'] - stats['zmin']
boxsize = max(xsize, ysize, zsize)

# half box size
boxhalf = boxsize / 2.0

# CHARMM scripting: crystal define cubic @boxsize @boxsize @boxsize 90 90 90
crystal.define_cubic(boxsize)
# CHARMM scripting: crystal build cutoff @boxhalf noper 0
crystal.build(boxhalf)

# Turn on image centering - bysegment for peptide, by residue for solvent
# CHARMM scripting: image byseg xcen 0 ycen 0 zcen 0 select segid adp end
image.setup_segment(0.0, 0.0, 0.0, 'ADP')
# CHARMM scripting: image byres xcen 0 ycen 0 zcen 0 select resname tip3 end
image.setup_residue(0.0, 0.0, 0.0, 'TIP3')

# Now specify nonbonded cutoffs for solvated box
cutnb = min(boxhalf,12)
cutim = cutnb
ctofnb = cutnb - 1.0
ctonnb = cutnb - 3.0

# Another nbonds example
# CHARMM scripting: nbonds cutnb @cutnb cutim @cutim ctofnb @ctofnb ctonnb @ctonnb -
#        inbfrq -1 imgfrq -1
pycharmm.NonBondedScript(
    cutnb=cutnb, cutim=cutim, ctonnb=ctonnb, ctofnb=ctofnb,
    eps=1.0,
    cdie=True,
    atom=True, vatom=True,
    fswitch=True, vfswitch=True,
    inbfrq=-1, imgfrq=-1).run()

# Fix the peptide and minimize the solvent to "fit"
# CHARMM scripting: cons fix select segid adp end
cons_fix.setup(pycharmm.SelectAtoms(seg_id='ADP'))

# Minimize the solvent positions with periodic boundary conditions using steepest descents
# CHARMM scripting: mini sd nstep 200 tole 1e-3 tolgrd 1e-3
minimize.run_sd(nstep=200, tolenr=1e-3, tolgrd=1e-3)

# Turn off fixed atoms
# CHARMM scripting: cons fix select none end
cons_fix.turn_off()

# Write the psf and coordinates for the solvated peptide
# write psf card name pdb/adp+wat.psf
write.psf_card('pdb/adp+wat.psf')
# write coor pdb name pdb/adp+wat_min.pdb
write.coor_pdb('pdb/adp+wat_min.pdb')