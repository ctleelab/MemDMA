#.gro to .pdb file conversion
# necessary for running original PackMem code as it takes only .pdb
from pathlib import Path

grofile = Path("/scratch/local/casakurai/packmem/extended-data/DPPC-DPPS-flat-production7.9-stripped.gro")
xtcfile = Path("/scratch/local/casakurai/packmem/extended-data/DPPC-DPPS-flat-production7.9-stripped.xtc")
PackMem_test_dir = Path("/scratch/local/casakurai/packmem/PackMem_original/testing-1800ext/DPPC-DPPS")
#PackMem_test_dir.mkdir()
interval = 5 #matching new PackMem frame analysis
#interval = 100 #2ns per frame
#interval = 10 #200ps per frame
def grotopdb(grofile, xtcfile):
    import MDAnalysis as mda
    u = mda.Universe(str(grofile),str(xtcfile))

    #want to write out a pdb file of each frame in the .xtc file
    for ts in u.trajectory[::interval]:
        u.atoms.write(f"{PackMem_test_dir}/pdbfile{ts.frame:05d}.pdb")

grotopdb(grofile,xtcfile)


