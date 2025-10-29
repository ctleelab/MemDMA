from pathlib import Path

list_types_atoms = ["W", "NA", "CL", "NC3", "PO4", "GL1", "GL2", "C1A", "D2A", "C3A", "C4A","C1B", "D2B", "C3B", "C4B", "C2A", "C2B" ]
list_types_atoms_lipid = ["NC3", "PO4", "GL1", "GL2", "C1A", "D2A", "C3A", "C4A","C1B", "D2B", "C3B", "C4B", "C2A", "C2B" ]
SIZE= .5
lipid = "name NC3 PO4 GL1 GL2 C1A D2A C3A C4A C1B D2B C3B C4B C2A C2B"
membrane_sel = "resname DPPC DOPC"


trajectory_path = Path("/home/casakurai/scratch/packmem/simulation_data/extracted_frames/DPPC_strain.2/defect_analysis")
