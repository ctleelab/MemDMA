import MDAnalysis as mda 
import subprocess 
from pathlib import Path

system = "DPPC_strain.2"
trajectory_path = Path(f"/home/casakurai/scratch/packmem/simulation_data_centered/{system}")


def _strip_trajectory(gro, xtc) -> None:

    # Input .gro (already provided)

    # Load your system
    u = mda.Universe(gro, xtc)

    # Select only lipid residues (replace with your lipid names)
    lipids = u.select_atoms("resname DPPC POPC DOPC")  

    # Write to a new .gro file
    lipids.write(f"{trajectory_path}/{system}_centered_lipid_only.gro")


    # Write a stripped trajectory (all frames)
    with mda.Writer(f"{trajectory_path}/{system}_centered_lipid_only.xtc", lipids.n_atoms) as W:
        for ts in u.trajectory:
            W.write(lipids)


_strip_trajectory(f"{trajectory_path}/production7.7.gro", f"{trajectory_path}/production7.7_centered_highest_lipid.xtc")