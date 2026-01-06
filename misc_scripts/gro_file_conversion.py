#.gro to .pdb file conversion
from pathlib import Path

grofile = Path("/home/casakurai/scratch/packmem/tests/DPPC_flat_lipid_only.gro")

def grotopdb(grofile):
    import MDAnalysis as mda
    u = mda.Universe(grofile)
    u.atoms.write("pdbfile.pdb")

    return "pdbfile.pdb"

grotopdb(grofile)


