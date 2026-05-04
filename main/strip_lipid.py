#!/usr/bin/env python3
import sys
import MDAnalysis as mda 
from pathlib import Path
import util

# pixi inputs
gro_file = sys.argv[1]
xtc_file = sys.argv[2]
file_name = str(sys.argv[3])
output_path = Path(sys.argv[4])

def strip_trajectory(gro,xtc,file_name, output_path) -> None:

    output_file_name = output_path / f"{file_name}_stripped"

    # Load your system
    u = mda.Universe(gro,xtc)

    # Select only lipid residues (replace with your lipid names)
    lipids = u.select_atoms(util.lipids)  

    # Write to a new .gro file
    lipids.write(f"{output_file_name}.gro")

    # Write a stripped trajectory (coordinates for only lipids)
    with mda.Writer(f"{output_file_name}.xtc", lipids.n_atoms) as W:
        for ts in u.trajectory:
            W.write(lipids)

    

strip_trajectory(gro_file, xtc_file,file_name,output_path)