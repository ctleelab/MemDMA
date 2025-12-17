#!/bin/bash

# Script to launch PackMem analysis on multiple PDB files.
# P. Fuchs, R. Gautier, 2016

# set prefix of the pdb files
prefix=""

# set PackMem absolute path
PackMemPATH=/Users/ctlee/PackMem-leelab/main
trajectoryPATH=/Users/ctlee/PackMem-leelab/simulation_data/

# # For looping over frames
# # set the total number of frames
# total_number_of_frames=0
# for pdbnum in $(seq 0 ${total_number_of_frames})


#set certain frame to be evaluated
pdbnum=200

    # print counter to screen
    echo "$(date): PackMem running on ${trajectoryPATH}/${pdbnum}.gro"
    # launch PackMem for the 3 types of packing defects
    ${PackMemPATH}/base_clean_v2.py -i ${trajectoryPATH}/${pdbnum}.gro \
                              -r ${PackMemPATH}/vdw_radii_CGMartini3_openbeta2_inaccuratedata_testing.txt \
                              -p ${PackMemPATH}/param_CGMartini3_openbeta2.txt \
                              -o ${pdbnum} -d 1.0 -t deep \
                              -f ${pdbnum}
