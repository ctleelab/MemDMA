from pathlib import Path
import MDAnalysis as mda

lipid = "DOPC-buckled"
trajectory_path = Path(
    f"/scratch/local/casakurai/packmem/tests"
)


def extract(gro_file, trajectory_file,stop_frame,stride):
    '''
    Creates a new .gro and .xtc file where the trajectory file is only the first 100 frames.
    
    Args:
    gro_file: A .gro file
    trajectory_file: A .xtc file
    stop_frame: number of frames to be included in new .xtc
    stride: The interval of frames written.

    Returns:
    lipid_trimmed.xtc: An .xtc file with stop_number of frames read at interval of stride

    '''


    u = mda.Universe(gro_file, trajectory_file)

    with mda.Writer(f"{lipid}_trimmed.xtc", n_atoms=u.atoms.n_atoms) as W:
        for ts in u.trajectory[:stop_frame:stride]:   # first 1000 frames
            W.write(u.atoms)


extract(
    trajectory_path / f"{lipid}.gro",
    trajectory_path / f"{lipid}.xtc",
    100,
    5
)