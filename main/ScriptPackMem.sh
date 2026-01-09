#!/bin/bash

# set prefix of the pdb files
prefix=""

PackMemPATH=.
trajectoryPATH=../tests

#set certain frame to be evaluated
#pdbnum=200
shape=DPPC_strain2_all_atoms

    # print counter to screen
    echo "$(date): PackMem running on ${trajectoryPATH}/${shape}.gro"
    # launch PackMem for the 3 types of packing defects
    ${PackMemPATH}/base_clean_v2.py -i ${trajectoryPATH}/${shape}.gro \
                              -r vdw_radii_CGMartini3_openbeta2_updated.txt \
                              -o ${shape} -f ${shape} -d 1
