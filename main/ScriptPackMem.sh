#!/bin/bash

# Usage: bash run_packmem.sh <system_NAME>

# Prompt user for system name
read -p "Enter system name (e.g., DPPC_strain.2): " system

# Check input
if [ -z "$system" ]; then
    echo "Error: No system name provided."
    exit 1
fi

export PACKMEM_SYSTEM="$system"

PackMemPATH=/home/casakurai/scratch/packmem/main
trajectoryPATH=/home/casakurai/scratch/packmem/simulation_data_centered/${system}/extracted_gro_frames


# set prefix of the pdb files
prefix=""

# For looping over frames
# set the total number of frames
total_number_of_frames=30000
for framenum in $(seq 0 1000 ${total_number_of_frames}) #increment by 10 frames each time
do

    # print counter to screen
    echo "$(date): PackMem running on ${trajectoryPATH}/${framenum}.gro"
    # launch PackMem for the 3 types of packing defects
    ${PackMemPATH}/base_clean_v2.py -i ${trajectoryPATH}/${framenum}.gro \
                                -r ${PackMemPATH}/vdw_radii_CGMartini3_openbeta2_updated.txt \
                                -p ${PackMemPATH}/param_CGMartini3_openbeta2.txt \
                                -o ${framenum} -d 0.1 \
                                -f ${framenum}

done
