#extract gro frames from .xtc file 

from pathlib import Path
trajectory_path = Path("/scratch/casakurai/packmem/simulation_data/extracted_frames/DPPC_strain.2")
files = Path("/scratch/casakurai/packmem/simulation_data")

def extract(topology_file, trajectory_file):
    import MDAnalysis as mda

    u = mda.Universe(topology_file, trajectory_file)

    # Save every frame as a .gro file
    # takes every 10th frame from the trajectory
    for i, ts in enumerate(u.trajectory[::10]):  # change stride as needed
        u.atoms.write(trajectory_path/f"{i}.gro")


extract(files/f"DPPC_strain2_equilibration7.6.gro", files/f"DPPC_strain2_equilibration7.6.xtc")