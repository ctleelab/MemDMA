#!/bin/bash

# Script to launch PackMem analysis on multiple PDB files.
# P. Fuchs, R. Gautier, 2016

# set prefix of the pdb files
prefix="pdbfile"

# set PackMem absolute path
# change this path according to your system 
#PackMemPATH=/home/clifford/fuchs/test_PackMem/packmem
read -p "Enter lipid composition: " lipid

PackMemPATH=/scratch/local/casakurai/packmem/PackMem_original
trajectoryPATH=/scratch/local/casakurai/packmem/PackMem_original/testing-1800ext/${lipid}
analysisPATH=/scratch/local/casakurai/packmem/PackMem_original/testing-1800ext/${lipid}/analysis/seperated-leaflets



# set the total number of frames
# in this example we have membrane0.pdb, membrane1.pdb, ..., membrane400.pdb

###################
# No change below #
###################

# The following Total_*.txt files will contain the statistics for each type of
# packing defects accumulated over all frames (one file per membrane leaflet)
# Delete those files if they already exist before launching the main loop
# (the -f flag avoids an error if the file doesn't exist)
rm -f Total_Up_${prefix}_deep.txt
rm -f Total_Lo_${prefix}_deep.txt
rm -f Total_Up_${prefix}_shallow.txt
rm -f Total_Lo_${prefix}_shallow.txt
rm -f Total_Up_${prefix}_all.txt
rm -f Total_Lo_${prefix}_all.txt

# loop over all frames
# for example, if we have 400ns with one frame every 1000 ps,
# we loop from 0 to 400 frames (0-400 ns)
for pdbfile in ${trajectoryPATH}/pdbfile*.pdb; do
    base=$(basename "${pdbfile}" .pdb)

    pdbnum=${base#pdbfile}

    echo "PDB number: ${pdbnum}"

    echo "$(date): PackMem running on ${pdbfile}"

    ${PackMemPATH}/PackMem-original.py \
        -i "${pdbfile}" \
        -r "${PackMemPATH}/vdw_radii_CGMartini3_openbeta2_updated.txt" \
        -p "${PackMemPATH}/param_CGMartini3_openbeta2.txt" \
        -o "${analysisPATH}${base}" \
        -op "${analysisPATH}" \
        -num "${pdbnum}" \
        -d 1.0 \
        -t deep
    ${PackMemPATH}/PackMem-original.py \
        -i "${pdbfile}" \
        -r "${PackMemPATH}/vdw_radii_CGMartini3_openbeta2_updated.txt" \
        -p "${PackMemPATH}/param_CGMartini3_openbeta2.txt" \
        -o "${analysisPATH}${base}" \
        -op "${analysisPATH}" \
        -num "${pdbnum}" \
        -d 1.0 \
        -t shallow
    ${PackMemPATH}/PackMem-original.py \
        -i "${pdbfile}" \
        -r "${PackMemPATH}/vdw_radii_CGMartini3_openbeta2_updated.txt" \
        -p "${PackMemPATH}/param_CGMartini3_openbeta2.txt" \
        -o "${analysisPATH}${base}" \
        -op "${analysisPATH}" \
        -num "${pdbnum}" \
        -d 1.0 \
        -t all

    # accumulate packing defects of the current frame in Total_*.txt files
    cat  ${analysisPATH}${base}_Up_Deep_result.txt >> ${analysisPATH}Total_Up_deep.txt
    cat  ${analysisPATH}${base}_Lo_Deep_result.txt >> ${analysisPATH}Total_Lo_deep.txt
    cat  ${analysisPATH}${base}_Up_Shallow_result.txt >> ${analysisPATH}Total_Up_shallow.txt
    cat  ${analysisPATH}${base}_Lo_Shallow_result.txt >> ${analysisPATH}Total_Lo_shallow.txt
    cat  ${analysisPATH}${base}_Up_All_result.txt >> ${analysisPATH}Total_Up_all.txt
    cat  ${analysisPATH}${base}_Lo_All_result.txt >> ${analysisPATH}Total_Lo_all.txt
    # we no longer need the defects of the current frame
    rm -f ${analysisPATH}${base}_Up_Deep_result.txt
    rm -f ${analysisPATH}${base}_Lo_Deep_result.txt
    rm -f ${analysisPATH}${base}_Up_Shallow_result.txt
    rm -f ${analysisPATH}${base}_Lo_Shallow_result.txt
    rm -f ${analysisPATH}${base}_Up_All_result.txt
    rm -f ${analysisPATH}${base}_Lo_All_result.txt
done

# # Finally, concatenate the packing defects statistics of the lower and upper leaflets
# cat ${analysisPATH}Total_Up_deep.txt ${analysisPATH}Total_Lo_deep.txt > ${analysisPATH}Total_deep.txt
# cat ${analysisPATH}Total_Up_shallow.txt ${analysisPATH}Total_Lo_shallow.txt > ${analysisPATH}Total_shallow.txt
# cat ${analysisPATH}Total_Up_all.txt ${analysisPATH}Total_Lo_all.txt > ${analysisPATH}Total_all.txt

# # Cleanup the statistics per leaflet
# rm -f ${analysisPATH}Total_Up_deep.txt ${analysisPATH}Total_Lo_deep.txt
# rm -f ${analysisPATH}Total_Up_shallow.txt ${analysisPATH}Total_Lo_shallow.txt
# rm -f ${analysisPATH}Total_Up_all.txt ${analysisPATH}Total_Lo_all.txt
