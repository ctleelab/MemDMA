import util
import MDAnalysis as mda 
import subprocess 


def _strip_trajectory(gro) -> None:

    # Input .gro (already provided)

    # Load your system
    u = mda.Universe(gro)

    # Select only lipid residues (replace with your lipid names)
    lipids = u.select_atoms("resname DPPC POPC DOPC")  

    # Write to a new .gro file
    lipids.write("lipid_only.gro")


_strip_trajectory("testing.gro")