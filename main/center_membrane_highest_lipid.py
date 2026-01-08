import MDAnalysis as mda
from pathlib import Path
import numpy as np

system="DOPC_strain.2"

#path
data=Path("/home/casakurai/scratch/packmem/simulation_data_centered/")
system_path = Path(f"/home/casakurai/scratch/packmem/simulation_data_centered/{system}")

#files 
trajectory= f"{system_path}/production7.7.tpr"
xtc = f"{system_path}/production7.7.xtc"
output_xtc = f"{system_path}/production7.7_centered_highest_lipid.xtc"

#load in the universe 
u = mda.Universe(trajectory, xtc)

#select lipids 
lipids = u.select_atoms("resname DOPC DPPC DOPS DPPS DOPA DPPA")
all_atoms = u.select_atoms("all")

with mda.Writer(output_xtc, all_atoms.n_atoms) as W:
    for ts in u.trajectory:

        z_positions = lipids.positions[:,2] #pulls out the z positions of all lipids
        peak_z_value = np.max(z_positions) #finds the max z position out of all lipids

        
        #select lipids within a cuttoff of the highest lipid
        z_threshold = peak_z_value - 1 #1A below the peak

        #create a boolean mask fo lipids that are within the radius
        peak_lipids_mask = z_positions <=  z_threshold

        #define the lpids that are within the peak region
        peak_lipids = lipids[peak_lipids_mask]

        #calcuate the center of mass of the peak lipids
        com_peak = peak_lipids.center_of_mass()

        #shift the box center based on the com of the peak lipids
        box_center = np.array([ts.dimensions[0]/2, ts.dimensions[1]/2, ts.dimensions[2]/2])
        shift = box_center - com_peak
        all_atoms.positions = all_atoms.positions.copy() + shift

        #PBC
        all_atoms.wrap(compound="atoms")

        #write out the frame
        W.write(all_atoms)











