#!/usr/bin/env python3
import sys
import MDAnalysis as mda 
from pathlib import Path
import numpy as np
import util

# pixi inputs
gro_file = sys.argv[1]
xtc_file = sys.argv[2]
file_name = str(sys.argv[3])
output_path = Path(sys.argv[4])



def centered_membrane(gro,xtc,file_name, output_path):
    output_file_name = output_path / f"{file_name}_centered"
    u = mda.Universe(gro, xtc)
    #select lipids 
    lipids = u.select_atoms(util.lipids)
    all_atoms = u.select_atoms("all")


    # --- First frame for GRO ---
    ts = u.trajectory[0]  # first frame

    # Center first frame based on peak lipids
    z_positions = lipids.positions[:, 2]
    x_positions = lipids.positions[:, 0]
    peak_z_value = np.max(z_positions) #finds the max z position out of all lipids
    peak_z_value_idx = np.argmax(z_positions)
    peak_lipid_x = lipids.positions[peak_z_value_idx, 0]
    
    z_threshold = peak_z_value - 4 #4nm below the peak, want it to be ~thickness of the membrane
    x_threshold_min = peak_lipid_x - 1
    x_threshold_max = peak_lipid_x + 1
    
    #create two boolean mask for lipids that are within the z threshold and for the x threshold
    #choose lipids that are within 10A of the peak z position
    peak_lipids_mask_z = z_positions >= z_threshold 
    peak_lipid_mask_x =  (x_positions >= x_threshold_min) & (x_positions <= x_threshold_max)

    #combined mask 
    total_peak_lipids_mask = peak_lipids_mask_z & peak_lipid_mask_x

    #define the lipids that are within the specified peak regions
    peak_lipids = lipids[total_peak_lipids_mask]

    #calculate the center of mass of the peak lipids
    com_peak = peak_lipids.center_of_mass()

    #shift the box center based on the com of the peak lipids
    box_center = np.array([ts.dimensions[0]/2, ts.dimensions[1]/2, ts.dimensions[2]/2]) #center the com at (0,0,0)
    shift = box_center - com_peak
    all_atoms.positions = all_atoms.positions.copy() + shift

    #PBC
    all_atoms.wrap(compound="atoms")


    # --- Save first frame as GRO ---
    all_atoms.write(f"{output_file_name}.gro")

    
    with mda.Writer(f"{output_file_name}.xtc", all_atoms.n_atoms) as W:
        for ts in u.trajectory:
            z_positions = lipids.positions[:,2] #pulls out the z positions of all lipids
            x_positions = lipids.positions[:,0] #pulls out the x positions of all lipids
            peak_z_value = np.max(z_positions) #finds the max z position out of all lipids

            #finds the indx of the the lipid wiht the hightest z value
            peak_z_value_idx = np.argmax(z_positions)

            #pulls out the (x,y) of the max lipid
            peak_lipid_x = lipids.positions[peak_z_value_idx, 0]

            #select lipids within a cuttoff of the highest lipid's z-coord
            z_threshold = peak_z_value - 4 #4nm below the peak, want it to be ~thickness of the membrane

            #select lipids within a cutoff of the highest lipid's x-coord to avoid two peaks being selected
            # +/- 1nm in the x direction
            x_threshold_min = peak_lipid_x - 1
            x_threshold_max = peak_lipid_x + 1

            #create two boolean mask for lipids that are within the z threshold and for the x threshold
            #choose lipids that are within 10A of the peak z position
            peak_lipids_mask_z = z_positions >= z_threshold 
            peak_lipid_mask_x =  (x_positions >= x_threshold_min) & (x_positions <= x_threshold_max)

            #combined mask 
            total_peak_lipids_mask = peak_lipids_mask_z & peak_lipid_mask_x

            #define the lipids that are within the specified peak regions
            peak_lipids = lipids[total_peak_lipids_mask]

            #calculate the center of mass of the peak lipids
            com_peak = peak_lipids.center_of_mass()

            #shift the box center based on the com of the peak lipids
            box_center = np.array([ts.dimensions[0]/2, ts.dimensions[1]/2, ts.dimensions[2]/2]) #center the com at (0,0,0)
            shift = box_center - com_peak
            all_atoms.positions = all_atoms.positions.copy() + shift

            #PBC
            all_atoms.wrap(compound="atoms")

            #write out the frame
            W.write(all_atoms)

    
centered_membrane(gro_file, xtc_file,file_name,output_path)