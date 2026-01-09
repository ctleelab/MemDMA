from pathlib import Path
import os 


system = os.environ.get("PACKMEM_SYSTEM") #gets the system name from environment variable


if system is None:
    raise RuntimeError(
        "PACKMEM_SYSTEM environment variable not set.\n"
        "Example: export PACKMEM_SYSTEM=DPPC_strain.2"
    )

list_types_atoms = [
    "W",
    "NA",
    "CL",
    "CNO"
    "NC3",
    "PO4",
    "GL1",
    "GL2",
    "C1A",
    "D2A",
    "C3A",
    "C4A",
    "C1B",
    "D2B",
    "C3B",
    "C4B",
    "C2A",
    "C2B",
]
list_types_atoms_lipid = [
    "CNO"
    "NC3",
    "PO4",
    "GL1",
    "GL2",
    "C1A",
    "D2A",
    "C3A",
    "C4A",
    "C1B",
    "D2B",
    "C3B",
    "C4B",
    "C2A",
    "C2B",
]
SIZE = 0.5
lipid = "name NC3 PO4 GL1 GL2 C1A D2A C3A C4A C1B D2B C3B C4B C2A C2B"
membrane_sel = "resname DPPC DOPC DOPA DOPS DPPA DPPS"

trajectory_path = Path(f"/home/casakurai/scratch/packmem/simulation_data_centered/{system}/extracted_gro_frames")
analysis_path = Path(f"/home/casakurai/scratch/packmem/simulation_data_centered/analysis_{system}/extracted_gro_frames")
