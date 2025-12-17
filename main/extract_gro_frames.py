# extract gro frames from .xtc file

from pathlib import Path

lipid = "DPPC"
trajectory_path = Path(
    f"/scratch/casakurai/packmem/simulation_data/extracted_frames/{lipid}_strain.2/prod"
)
files = Path("/scratch/casakurai/packmem/simulation_data")


def extract(topology_file, trajectory_file):
    import MDAnalysis as mda

    u = mda.Universe(topology_file, trajectory_file)

    # Save every frame as a .gro file
    # takes every 10th frame from the trajectory
    for i, ts in enumerate(u.trajectory[::25]):  # change stride as needed
        u.atoms.write(trajectory_path / f"{i}.gro")


extract(
    files / f"{lipid}_strain2_production7.7.gro",
    files / f"{lipid}_strain2_production7.7.xtc",
)
