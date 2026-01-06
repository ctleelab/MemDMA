import MDAnalysis as mda
import pdb
from pathlib import Path 

path = Path("/home/casakurai/scratch/packmem/tests")
u = mda.Universe(f"{path}/DOPC_flat_lipid_only.gro")

# pdb.set_trace()

new_resid_map = {
    old_id: new_id for new_id, old_id in enumerate(u.residues.resids, start=1)
}


for res in u.residues:
    res.resid = new_resid_map[res.resid]

with mda.Writer("out.gro", reindex=True) as w:
    w.write(u.atoms)

with mda.Writer("out.pdb") as w:
    w.write(u.atoms) 