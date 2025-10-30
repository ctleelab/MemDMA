#!/usr/bin/env python3

#general libraries
import argparse
import sys
import pdb as pdb
import numpy as np
import re 
import MDAnalysis as mda 
import pickle
from scipy.spatial import cKDTree
import csv


#library for showing progress bar for runs
from tqdm.auto import tqdm
#library that does parallel processing 
from tqdm.contrib.concurrent import process_map 


#importing my personal scripts
import param as p
import BasicFunctions as bfrg
import dico as d
import listes as l
import leaflet_finder as lf
import quadrilateral_mesh as qm
import mesh_clean as m
import util as util
import grid_slicing as gs
import preprocessing as pp 
import ray_tracing_parallelization_v3 as ry
import quad_centroid as qc
import neighborhood_grids as ng
import ray_tracing_quads as rq
import atom_grid_mapping as agm
import neighborhood_grid_graph as neigh_graph


##this chunk of code will only be run if the script is executed directly
##this is the interface with the user, so that they can specify inputs
if __name__ == "__main__":
    
    try:
        outputname = "output"
        pdbout = 0
        ##this is telling the flags of the program
        ##-i: tells you the flag that the user will type
        ##action: tells python what do the flag input given by the user
        ##dest: stores the file given by the user in the varible "filename", in args. So when you type arg.filename you get the file give by the user
        parser = argparse.ArgumentParser(description = 'Arguments for the app')
        ##changed to take .gro files
        parser.add_argument('-i', action = 'store', dest = 'filename', help = 'GRO File')
        parser.add_argument('-r', action = 'store', dest = 'filesrad',
                            help = 'File for radius (default vdw_radiiFinal2014.txt)',
                            default = 'vdw_radiiFinal2014.txt')
        parser.add_argument('-o', action = 'store', dest = 'outputname',
                            default = 'output', help = 'Name for output file (default output)')
        parser.add_argument('-d', action = 'store', dest = 'dist_suppl_Z', default = 1.0, type = float,
                            help = '(default value 1.0)')
        parser.add_argument('-t', action = 'store', dest = 'pd_type', default = 'all',
                            help = 'Packing Defect Type (all/shallow/deep) default all')
        parser.add_argument('-p', action = 'store', dest = 'paramFile', default = 'param.txt',
                            help = 'File for lipid parameters')
        parser.add_argument('-n', action = 'store', dest = 'indexFile',
                            help = 'File for index (Gromacs ndx style). Only Lower/Upper group accepted')
        parser.add_argument('-v', dest = 'pdbout', action = 'store_true',
                            help = 'Increase the verbosity')
        parser.add_argument('-f', action = 'store', dest = 'frame_num',
                            help = 'frame number')
        
        ##parse_args() will read the arguments that the user gives. These agruments were added to parser. This line of code then stores this information into the python variable args.        
        args = parser.parse_args()
        
        ##This is just front-end error handling with the user. Ensure that the user gives a pdb file
        if len(sys.argv) == 1:  ##sys.argv is the list of all arguments that the user gives when executing the script
            parser.print_usage() ##tells the user what arguments are expected
            sys.exit()

        ##takes in the .txt file (which describes the name of the lipid head beads and GL2)
        ##will output 3 dictionaries containing 
        # 1) DICT_3L taking lipid_bead head name: lipid name(3-letters) ex:"{'DPPC_NC3': 'DPPC', 'DOPC_NC3': 'DOPC'}"
        # 2) RESNAME_GLYC lipid name: GL2 bead head ex:"{'DPPC': 'GL2', 'DOPC': 'GL2'}"
        # 3) LIPID lipid name (3-letters) ex:"['DPPC', 'DOPC']"
        DICT_3L, RESNAME_GLYC, LIPID = p.set_params(args.paramFile)



        # determine packing defect type (all/shallow/deep)
        if args.pd_type == 'all':
            FlagPDtype = 0
        elif args.pd_type == 'shallow':
            FlagPDtype = 2
        elif args.pd_type == 'deep' :
            FlagPDtype = 1
        else: 
            print('ERROR : Packing defect type not known. Please choose : all/deep/shallow')
            sys.exit()
    
        ##error handling if the choosen distance is negative. This distance is key when determining where each bead falls within the mesh grid. 
        if args.dist_suppl_Z < 0.0 :
            print('ERROR : The distance for the -d option must be > 0')
            sys.exit()
        
        grolines = bfrg.read_file(args.filename) ##reads line of the pdb file and stores each line  
        radius = bfrg.read_radius(args.filesrad) ##creates a dictionary that stores lipid_atom:radius of atom
        
        aliphatic = bfrg.read_aliphatic(args.filesrad) ##creates a dictionary that stores lipid_atom: n vs a (polar vs aliphatic)
    except:
        print('Command line: PackMem.py -i file.pdb -r fileRadius.txt \
            -p param.txt -o output -d distGlyc -t deep/all/shallow [-n index file] [-v] or -h for help')
        sys.exit()

    #extract num_frame in pdb file (line MODEL) or 1 by default
    ##one per each pdb file, because you feed in pdb files one at a time 
    num_frame=pdb.find_numframe(grolines)


    #####################################################################
    ########## create a quadrilateral mesh of the membrane ##############
    #####################################################################

    #inputs: grofile, selection of beads for mesh, step size for mesh grid (nm)
    #outputs: (X, Y, Z_upper, Z_lower) these are all coordinates 
    ## leaflet_residues_new: {upper:[1,2,3], lower:[4,5,6]}
    ##note this also calculates the mesh information

    #Is the mesh is in nm?
    #Right now the mesh is 5A
    print("Start: running mesh")
    X, Y, zupper, zlower, leaflet, leaflet_residues, z_midplane = qm.mesh(args.filename, "name GL2", .1)
    print("Finished:running mesh")

    
    ##formating data for analysis
    print("Start:formatting data")
    X_f, Y_f, zupper_f, zlower_f = m.format_data(X,Y, zupper,zlower)
    print("Finished:formating data")


    ####################################################
    ########## Determining atom features ###############
    ####################################################

    print("Start:determining atom features")
    ##input: grolines, list of the types of atoms in the gro file, which residues are in the upper vs lower leaflet
    ##output: {atm#: GL2', 'lower', array([18.222,  4.37 , 14.85 ]), 'DOPC GL2', resnumber }

    list_types_atoms = util.list_types_atoms_lipid  #just loads in the types of atoms in the .gro file
    atom_features = pp.atom_features(grolines, list_types_atoms,leaflet_residues) 

    print("Finished:determining atom features")

    # #################################################################################
    # ########## compute the normal vector of each triangle in the mesh ###############
    # #################################################################################

    #inputs: X_f, Y_f, zupper_f, zlower_f
    #outputs: For the upper and lower leaflet of the membrane: centroids of the triangles, normal vectors of each triangle, 3 points that create each triangle, and the edges that compose each triangle
    #centroid_triangles_u {grid_number: corods of centroidxs (x,y,z)}
    ##note this also calculates the mesh information
    

    print("Start:computing normal vectors of triangle")
    centroid_triangles_u, centroid_triangles_l, norm_vector_triangles_up, norm_vector_triangles_l, coords_triangles_upper, coords_triangles_lower, edges_triangles_upper, edges_triangles_lower, quadrilateral_normal_upper, quadrilateral_normal_lower, centroid_triangles_lower_matrix, centroid_triangles_upper_matrix = m.normal_vector_triangles(X_f,Y_f,zupper_f,zlower_f)
    print("Finished:computing normal vectors of triangle")

   


    neigh_graph.plot_triangle_grid_assignment(centroid_triangles_u, centroid_triangles_l)

    #################################################################################
    ########## Determine the grid closest to the GL2 bead of each lipid #############
    #################################################################################

    #outline of code
    #1) pull out gl2 bead and its coordinates
    #2) create a neighborhood by cKDtree
    #3) determine which grid is closest to the GL2 bead by euclidean distance
    #5) assign a neighborhood to each bead based on lipid it belongs too. 

    print("Start:GL2 grid neighborhood assignment")
    #1-1) creating dictionaries for gl2 beads and all atoms
    #  gl2_beads_upper a dictionary of just gl2 beads {atm#: coords}
    #  gl2_beads_upper_resnumber a dictionary of just gl2 beads {atm#: resnumber}
    gl2_beads_upper = {}
    gl2_beads_lower = {}
    gl2_beads_upper_resnumber = {}
    gl2_beads_lower_resnumber = {}
    for key, v in atom_features.items(): #need .items if there are multiple items in the value of the key of a dictionary 
        if v[0] == "GL2":
            if v[1] == "upper":
                gl2_beads_upper[key] = v[2]*10 #add coords of gl2 beads but multiply by 10 to bring units to A
                gl2_beads_upper_resnumber[key] = v[4]
            if v[1] == "lower":
                gl2_beads_lower[key] = v[2]*10 #add coords of gl2 beads but multiply by 10 to bring units to A
                gl2_beads_lower_resnumber[key] = v[4]


    
    ##1-2) all atom dictionary
    ## output {atom number: res_number, coords} {3183': (12053, array([ 72.14,  58.45, 194.69]))}
    all_atoms_upper = {}
    all_atoms_lower = {}
    for key, v in atom_features.items(): #need .items if there are multiple items in the value of the key of a dictionary 
        if v[1] == "upper":
            all_atoms_upper[key] = (v[4], v[2]*10) #adds resnumber and coords(A, that's why multiply by 10)
        if v[1] == "lower":
            all_atoms_lower[key] = (v[4], v[2]*10) #adds resnumber and coords(A, that's why multiply by 10)

    #1-2) turn dictionary into a matrix, and have the keys. The index matches the keys to the matrix
    #gl2_beads_upper_keys = {k: None for k in gl2_beads_upper.keys()}
    #gl2_beads_upper_matrix = np.array(list(gl2_beads_upper.values()))
    #gl2_beads_lower_keys = {k: None for k in gl2_beads_lower.keys()}
    #gl2_beads_lower_matrix = np.array(list(gl2_beads_lower.values()))


    ##########################################
    ##calculations for quadrilateral grids####
    ##########################################
    # #2-1) create a matrix of the center of mass of each mesh grid for quadrilaterals
    #upper_centroids, lower_centroids, centroids_upper_dict, centroids_lower_dict  = qc.centroid_quadrilaterals(X_f, Y_f, zupper_f, zlower_f)
    # #2-2a) ##create neighborhood for each grid that has an upper limit of 10A with the quadrilaterals
    # upper_neighborhood = ng.neighborhood(centroids_upper_dict, 1)
    # lower_neighborhood = ng.neighborhood(centroids_lower_dict, 1)


    # #2) ##create neighborhood for each grid that has an upper limit of 10A with the triangles
    # ##This takes alot of memory 
    # print("Start: neighborbood")
    # upper_neighborhood = ng.neighborhood(centroid_triangles_u, 10)

    # # #get more information about the neighborhood created
    # # with open ("upper_neighborhood_cKDtree.csv", "w", newline= "") as f:
    # #     writer = csv.writer(f)
    # #     for key, value in upper_neighborhood.items():
    # #         writer.writerow([key, value])

    # lower_neighborhood = ng.neighborhood(centroid_triangles_l, 10)
    # print("Finished: neighborbood")

    # ##############################################
    # ######troubleshooting neighborhood grids######
    # ##############################################
    # if 50 in upper_neighborhood:
    #     print(upper_neighborhood[50])
    # else:
    #     print("not found")

    # if 100 in upper_neighborhood:
    #     print(upper_neighborhood[100])
    # else:
    #     print("not found")

    #######################################
    ###create a cKDtree of the centroids###
    #######################################
    tree_upper, grid_keys_upper, centroids_upper = ng.cKDtree_grid_neighborhood(centroid_triangles_u)

    tree_lower, grid_keys_lower, centroids_lower = ng.cKDtree_grid_neighborhood(centroid_triangles_l)


    #3) calculate the euclidean distance between the Gl2 and the grids in the mesh
    # can compare the array of Gl2 coordinates to the array of mesh COM (this should be quick because can do matrix math). Will compare by euclidean distance. 

    #inputs will be the centroids of each grid and the gl2 bead coordinates
    #output : {'4': (np.int64(42), array([ 8.09204793,  1.12890854, 19.91730258]}
    #output : {gl2 atom_number: grid_number, (coords of grid), list of neighborhood grids}
    #for some reason the gl2_beads_upper is in nm 

    #rather than computing the entire neighborhood for each grid, only compute the neighborhood for those that are assigned 

    print("Start: gl2 grid assignment")
    gl2_grid_assignment_upper = qc.grid_assignment(gl2_beads_upper, centroid_triangles_upper_matrix, centroid_triangles_u, tree_upper, grid_keys_upper, centroids_upper, 10) 
    gl2_grid_assignment_lower = qc.grid_assignment(gl2_beads_lower, centroid_triangles_lower_matrix, centroid_triangles_l, tree_lower, grid_keys_lower, centroids_lower, 10)
    print("Finished: gl2 grid assignment")

    neigh_graph.plot_gl2_grid_assignment(gl2_grid_assignment_upper, gl2_grid_assignment_lower,gl2_beads_lower)

    # #############################################################
    # ########## assign the neighborhood to each atom  ############
    # #############################################################

    #note: this is taking alot of memory...
    # #4) 
    #goal: assign to each atom a neighborhood of grids
    #matchs by res_number by using the Gl2 bead
    #output {atom number: res_number, coords, list of close by grids}

    atom_grid_neighborhood_assignment_upper = qc.atom_grid_mapping(all_atoms_upper, gl2_grid_assignment_upper, gl2_beads_upper_resnumber)
    atom_grid_neighborhood_assignment_lower = qc.atom_grid_mapping(all_atoms_lower, gl2_grid_assignment_lower, gl2_beads_lower_resnumber)

    print("Finished:GL2 grid neighborhood assignment")


    # # ##################################################################################################################
    # # ########## Determine the coordinate of an atom if projected onto the same plane as a given quadrilateral #########
    # # ##################################################################################################################
    
    # # #need the inputs of the neighborhood of grids 
    # # #inputs: quadrilateral_upper (which has the average normal vector and the quadrilateral number), atom_coords, upper_centroids (quadrilateral_centroids)
    # # #quadrilateral_upper {quad_number: coords of normal vector}
    # # #atom_grid_assignment_upper  {atom number: res_number, coords, list of close by grids}
    # # #upper_centroids is a matrix
    # # output = rq.project_block_of_triangles(quadrilateral_normal_upper, atom_grid_assignment_upper , upper_centroids)


    # # ##################################################################################################################
    # # ########## Determine the coordinate of an atom if projected onto the same plane as a given triangle      #########
    # # ##################################################################################################################
    
    print("Start: atom ray projection")

    ## grid_coords_np: is the (x,y,z) coordinates of the atoms in the upper leaflet
    ## grid_upper_num: is a list of grid numbers in the upper leaflet 
    grid_upper_num, grid_coords_upper_np = ry.extract_grid_coords(centroid_triangles_u)
    grid_lower_num, grid_coords_lower_np = ry.extract_grid_coords(centroid_triangles_l)

    atom_upper_num, atoms_coords_upper_np = ry.extract_atoms(atom_features, "upper")
    atom_lower_num, atoms_coords_lower_np = ry.extract_atoms(atom_features, "lower")
 

    #atom_grid_neighborhood_assignment_lower this needs to be in the input somewhere
    #output: neighborhood_grid_id: {atom number: list of grid number}
    #ouput: neighborhood_projections: {atom number: list of array(x,y,z)}
    #neighborhood_projections_upper, neighborhood_grid_id_upper = ry.neighborhood_project_atoms_to_triangles(atom_upper_num, grid_upper_num, grid_coords_upper_np, norm_vector_triangles_up, atom_grid_neighborhood_assignment_upper)
    #neighborhood_projections_lower, neighborhood_grid_id_lower = ry.neighborhood_project_atoms_to_triangles(atom_lower_num, grid_lower_num, grid_coords_lower_np, norm_vector_triangles_l, atom_grid_neighborhood_assignment_lower)

    print("Finished: atom ray projection")

    ########################################################
    ##  Determine if the atoms lies inside in the triangle##
    ########################################################
    
    print("Start: atom mapping to triangles")

    # with open("neighborhood_projections_upper_test.pickle", "rb") as f:
    #     neighborhood_projections_upper = pickle.load(f)

    # with open("neighborhood_projections_atom_id_test.pickle", "rb") as f:
    #     neighborhood_grid_id_upper = pickle.load(f)

    
    triangle_atom_assignment_upper  = {} 
    triangle_atom_assignment_lower  = {} 

    #trouble-shooting    
    #success = [] 
    #fails = [] 
    #print(max(neighborhood_grid_id_upper) )

    #output: neighborhood_grid_id: {atom number: list of grid number}
    #atom mapping to matrix quantifying polar vs nonpolar 
    #atom_grid_assignment_upper  {atom number: res_number, coords, list of close by grids}
    triangle_atom_assignment_upper = agm.sum_atom_mapping_upper(atom_features, radius, aliphatic, edges_triangles_upper, atom_grid_neighborhood_assignment_upper, coords_triangles_upper, grid_upper_num,  grid_coords_upper_np, norm_vector_triangles_up, args.frame_num )
    #triangle_atom_assignment_lower = agm.sum_atom_mapping_lower(atom_features, radius, aliphatic, edges_triangles_lower, neighborhood_grid_id_lower, coords_triangles_lower, grid_lower_num)

    print("completed frame")
    #Issue: too many grids are empty 

    #############################################
    #####plot the grid with the atom_mappings####
    #############################################

    #match grid number to coords of centroids
    #input: centroid_triangles_u and triangle_atom_assignment_upper 
    #I need x and y

    #first plot the x,y 

    #pull out the (x,y) from all triangle numbers


      
    