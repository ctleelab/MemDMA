#!/bin/bash

# Script to launch PackMem analysis on multiple PDB files.
# P. Fuchs, R. Gautier, 2016

# set prefix of the pdb files
prefix=""

# set PackMem absolute path
PackMemPATH=/home/casakurai/scratch/packmem/main
trajectoryPATH=/home/casakurai/scratch/packmem/tests


# # For looping over frames
# # set the total number of frames
# total_number_of_frames=0
# for pdbnum in $(seq 0 ${total_number_of_frames})


#set certain frame to be evaluated
#pdbnum=200
shape=DPPC_flat_lipid_only

    # print counter to screen
    echo "$(date): PackMem running on ${trajectoryPATH}/${shape}.gro"
    # launch PackMem for the 3 types of packing defects
    ${PackMemPATH}/base_clean_v2.py -i ${trajectoryPATH}/${shape}.gro \
                              -r ${PackMemPATH}/vdw_radii_CGMartini3_openbeta2_updated.txt \
                              -p ${PackMemPATH}/param_CGMartini3_openbeta2.txt \
                              -o ${shape} -d 0.1 \
                              -f ${shape}
