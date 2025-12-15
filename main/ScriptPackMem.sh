#!/bin/bash

# Script to launch PackMem analysis on multiple PDB files.
# P. Fuchs, R. Gautier, 2016

# set prefix of the pdb files
prefix=""

# set PackMem absolute path
PackMemPATH=/home/casakurai/scratch/packmem/main
trajectoryPATH=/home/casakurai/scratch/packmem/simulation_data/extracted_frames/DPPC_strain.2

# # set the total number of frames
# total_number_of_frames=0
# in this example we have membrane0.pdb, membrane1.pdb, ..., membrane400.pdb

###################
# No change below #
###################

# The following Total_*.txt files will contain the statistics for each type of
# packing defects accumulated over all frames (one file per membrane leaflet)
# Delete those files if they already exist before launching the main loop
# (the -f flag avoids an error if the file doesn't exist)
# rm -f Total_Up_${prefix}_deep.txt
# rm -f Total_Lo_${prefix}_deep.txt
# rm -f Total_Up_${prefix}_shallow.txt
# rm -f Total_Lo_${prefix}_shallow.txt
# rm -f Total_Up_${prefix}_all.txt
# rm -f Total_Lo_${prefix}_all.txt

# # loop over all frames
# # for example, if we have 400ns with one frame every 1000 ps,
# # we loop from 0 to 400 frames (0-400 ns)
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
    # ${PackMemPATH}/PackMem.py -i pdb/${prefix}${pdbnum}.pdb \
    #                           -r ${PackMemPATH}/vdw_radii_Charmm.txt \
    #                           -p ${PackMemPATH}/param_Charmm.txt \
    #                           -o ${prefix}${pdbnum} -d 1.0 -t shallow
    # ${PackMemPATH}/PackMem.py -i pdb/${prefix}${pdbnum}.pdb \
    #                           -r ${PackMemPATH}/vdw_radii_Charmm.txt \
    #                           -p ${PackMemPATH}/param_Charmm.txt \
    #                           -o ${prefix}${pdbnum} -d 1.0 -t all
    # accumulate packing defects of the current frame in Total_*.txt files
    # cat  ${prefix}${pdbnum}_Up_Deep_result.txt >> Total_Up_${prefix}_deep.txt
    # cat  ${prefix}${pdbnum}_Lo_Deep_result.txt >> Total_Lo_${prefix}_deep.txt
    # cat  ${prefix}${pdbnum}_Up_Shallow_result.txt >> Total_Up_${prefix}_shallow.txt
    # cat  ${prefix}${pdbnum}_Lo_Shallow_result.txt >> Total_Lo_${prefix}_shallow.txt
    # cat  ${prefix}${pdbnum}_Up_All_result.txt >> Total_Up_${prefix}_all.txt
    # cat  ${prefix}${pdbnum}_Lo_All_result.txt >> Total_Lo_${prefix}_all.txt
    # # we no longer need the defects of the current frame
    # rm -f ${prefix}${pdbnum}_Up_Deep_result.txt
    # rm -f ${prefix}${pdbnum}_Lo_Deep_result.txt
    # rm -f ${prefix}${pdbnum}_Up_Shallow_result.txt
    # rm -f ${prefix}${pdbnum}_Lo_Shallow_result.txt
    # rm -f ${prefix}${pdbnum}_Up_All_result.txt
    # rm -f ${prefix}${pdbnum}_Lo_All_result.txt


# # Finally, concatenate the packing defects statistics of the lower and upper leaflets
# cat Total_Up_${prefix}_deep.txt Total_Lo_${prefix}_deep.txt > Total_${prefix}_deep.txt
# cat Total_Up_${prefix}_shallow.txt Total_Lo_${prefix}_shallow.txt > Total_${prefix}_shallow.txt
# cat Total_Up_${prefix}_all.txt Total_Lo_${prefix}_all.txt > Total_${prefix}_all.txt

# # Cleanup the statistics per leaflet
# rm -f Total_Up_${prefix}_deep.txt Total_Lo_${prefix}_deep.txt
# rm -f Total_Up_${prefix}_shallow.txt Total_Lo_${prefix}_shallow.txt
# rm -f Total_Up_${prefix}_all.txt Total_Lo_${prefix}_all.txt
