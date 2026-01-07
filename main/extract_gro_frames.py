# extract gro frames from .xtc file

from pathlib import Path

system = "DPPC_strain.2"
trajectory_path = Path(f"/home/casakurai/scratch/packmem/simulation_data_centered/{system}")
extracted_frames_path = Path(f"/home/casakurai/scratch/packmem/simulation_data_centered/{system}/extracted_gro_frames")

files = Path("/scratch/casakurai/packmem/simulation_data")


def extract(topology_file, trajectory_file):
    import MDAnalysis as mda

    u = mda.Universe(topology_file, trajectory_file)

    # Save every frame as a .gro file
    # takes every 10th frame from the trajectory
    for i, ts in enumerate(u.trajectory[::25]):  # change stride as needed
        u.atoms.write(extracted_frames_path / f"{i}.gro")


extract(
    trajectory_path  / f"{system}_centered_lipid_only.gro",
    trajectory_path  / f"{system}_centered_lipid_only.xtc",
)
