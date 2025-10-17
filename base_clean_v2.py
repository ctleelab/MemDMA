#general libraries
import argparse
import sys
import pdb as pdb
import numpy as np
import re 
import MDAnalysis as mda 
import pickle
from scipy.spatial import cKDTree


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
import util
import grid_slicing as gs
import preprocessing as pp 
import ray_tracing_parallelization_v2 as ry
import matching 
import quad_centroid as qc
import neighborhood_grids as ng
import ray_tracing_quads as rq


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

    #mesh is in nm
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

    
    #################################################################################
    ########## Determine the grid closest to the GL2 bead of each lipid #############
    #################################################################################

    #outline of code
    #1) pull out gl2 bead and its coordinates
    #2) determine which grid is closest to the GL2 bead (can pull from information from CG2TS)
    #3) create a neighborhood around that grid
    #4) assign the grid to the lipid such that now every bead in that lipid gets assigned the grid
    
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
                gl2_beads_upper[key] = v[2]
                gl2_beads_upper_resnumber[key] = v[4]
            if v[1] == "lower":
                gl2_beads_lower[key] = v[2]
                gl2_beads_lower_resnumber[key] = v[4]


    
    ##all atom dictionary
    ## output {atom number: res_number, coords}
    all_atoms_upper = {}
    all_atoms_lower = {}
    for key, v in atom_features.items(): #need .items if there are multiple items in the value of the key of a dictionary 
        if v[1] == "upper":
            all_atoms_upper[key] = (v[4], v[2]) #adds resnumber 
        if v[1] == "lower":
            all_atoms_lower[key] = (v[4], v[2]) #adds resnumber 


    #1-2) turn dictionary into a matrix, and have the keys. The index matches the keys to the matrix
    #gl2_beads_upper_keys = {k: None for k in gl2_beads_upper.keys()}
    #gl2_beads_upper_matrix = np.array(list(gl2_beads_upper.values()))
    #gl2_beads_lower_keys = {k: None for k in gl2_beads_lower.keys()}
    #gl2_beads_lower_matrix = np.array(list(gl2_beads_lower.values()))



    #2-1) create a matrix of the center of mass of each mesh grid
    upper_centroids, lower_centroids, centroids_upper_dict, centroids_lower_dict  = qc.centroid_quadrilaterals(X_f, Y_f, zupper_f, zlower_f)


    #2-2) ##create neighborhood for each grid that has an upper limit of 10A 
    upper_neighborhood = ng.neighborhood(centroids_upper_dict, 1)
    lower_neighborhood = ng.neighborhood(centroids_lower_dict, 1)

    ##############################################
    ######troubleshooting neighborhood grids######
    ##############################################
    
    # if 1993802 in upper_neighborhood:
    #     print(upper_neighborhood[1993802])
    # else:
    #     print("not found")

    
    # if 100 in upper_neighborhood:
    #     print(upper_neighborhood[100])
    # else:
    #     print("not found")


    #2-3) calculate the euclidean distance between the Gl2 and the grids in the mesh
    # can compare the array of Gl2 coordinates to the array of mesh COM (this should be quick because can do matrix math). Will compare by euclidean distance. 

    #inputs will be the centroids of each grid and the gl2 bead coordinates
    #output : {'4': (np.int64(42), array([ 8.09204793,  1.12890854, 19.91730258]}
    #output : {gl2 atom_number: grid_number, (coords of grid), list of neighborhood grids}
    gl2_grid_assignment_upper = qc.grid_assignment(gl2_beads_upper, upper_centroids, upper_neighborhood) 
    gl2_grid_assignment_lower = qc.grid_assignment(gl2_beads_lower, lower_centroids, lower_neighborhood) 


   
  
    # #############################################################
    # ########## assign the neighborhood to each lipid ############
    # #############################################################

    # #3-1) 
    #input 
    #match by res_number 
    #output {atom number: res_number, coords, list of close by grids}

    atom_grid_assignment_upper = qc.atom_grid_mapping(all_atoms_upper, gl2_grid_assignment_upper, gl2_beads_upper_resnumber)
    atom_grid_assignment_lower = qc.atom_grid_mapping(all_atoms_lower, gl2_grid_assignment_lower, gl2_beads_lower_resnumber)

    print("Finished:GL2 grid neighborhood assignment")

    # #################################################################################
    # ########## compute the normal vector of each triangle in the mesh ###############
    # #################################################################################

    #inputs: X_f, Y_f, zupper_f, zlower_f
    #outputs: For the upper and lower leaflet of the membrane: centroids of the triangles, normal vectors of each triangle, 3 points that create each triangle, and the edges that compose each triangle
    ##note this also calculates the mesh information

    print("Start:computing normal vectors of triangle")
    centroid_triangles_u, centroid_triangles_l, norm_vector_triangles_up, norm_vector_triangles_l, coords_triangles_upper, coords_triangles_lower, edges_triangles_upper, edges_triangles_lower, quadrilateral_normal_upper, quadrilateral_normal_lower = m.normal_vector_triangles(X_f,Y_f,zupper_f,zlower_f) #potentially need to increase the efficiency of this code
    print("Finished:computing normal vectors of triangle")


    # ##################################################################################################################
    # ########## Determine the coordinate of an atom if projected onto the same plane as a given quadrilateral #########
    # ##################################################################################################################
    
    #need the inputs of the neighborhood of grids 
    #inputs: quadrilateral_upper (which has the average normal vector and the quadrilateral number), atom_coords, upper_centroids (quadrilateral_centroids)
    #quadrilateral_upper {quad_number: coords of normal vector}
    #atom_grid_assignment_upper  {atom number: res_number, coords, list of close by grids}
    #upper_centroids is a matrix
    output = rq.project_block_of_triangles(quadrilateral_normal_upper, atom_grid_assignment_upper , upper_centroids)





    # #calculate the tree = cKDTree(atom_coords) prior 
    # tree, atom_num, atom_coords_np = ry.tree_cKDTree(atom_features, "upper")


    # # print("Start: atom ray projection")
    # neighborhood_projections_upper, neighborhood_atom_id_upper = ry.neighborhood_project_atoms_to_triangles_tree(atom_num, atom_coords_np, centroid_triangles_u, norm_vector_triangles_up, tree, 6 )

    # with open("neighborhood_projections_upper_test.pickle", "wb") as f:
    #     pickle.dump(neighborhood_projections_upper , f) 

    # with open("neighborhood_projections_atom_id_test.pickle", "wb") as f:
    #     pickle.dump(neighborhood_atom_id_upper, f) 
    # print("Finished: atom ray projection")



    # #######################################################
    # #now determine what atoms lie inside in each triangle##
    # #######################################################
    
    # with open("neighborhood_projections_upper_v1.pickle", "rb") as f:
    #     neighborhood_projections_upper = pickle.load(f)

    # with open("neighborhood_projections_atom_id_v1.pickle", "rb") as f:
    #     neighborhood_atom_id_upper = pickle.load(f)

    
    # triangle_atom_assignment_upper  = {} 
    # triangle_atom_assignment_lower  = {} 

    # print("Start: atom mapping to triangles")

    # for triangle_num, atom_ids in enumerate(neighborhood_atom_id_upper): #loops over each triangle 
    #     stop_triangle = False
    #     #if triangle_num < 20000: #just emplaced so that I can troubleshoot the code
    #     print(triangle_num) #here to help troubleshoot the code 
    #     # Determine the coordinates for the atom by taking atom_id and searching in the same triangle
    #     # in neighborhood_projections and pulling out its coordinates
    #     for atom_id in atom_ids: #loops over the atoms in each triangle
    #         if stop_triangle:   #will move onto a new triangle if an polar atom was mapped
    #             break
    #         if atom_features[atom_id][1] == "lower": ##skips atoms that belong to the lower leaflet
    #             continue
    #         #coordinates of atom in plane of triangle 
    #         atom_coords = matching.get_projection(triangle_num, atom_id, neighborhood_projections_upper, neighborhood_atom_id_upper)

    #         #determine the radius of the atom you're working with
    #         lipid_atm_name = atom_features[atom_id][3]
    #         atm_radius = radius[lipid_atm_name]
    #         atm_type = aliphatic[lipid_atm_name]


    #         triangle_idex_num = int(triangle_num)+1

    #         ##import edges of triangle and import points of tringles
    #         edges_and_points = {
    #             "edge1": (edges_triangles_upper[triangle_idex_num][0], coords_triangles_upper[triangle_idex_num][0]),
    #             "edge2": (edges_triangles_upper[triangle_idex_num][1], coords_triangles_upper[triangle_idex_num][1]),
    #             "edge3": (edges_triangles_upper[triangle_idex_num][2], coords_triangles_upper[triangle_idex_num][2])
    #         }

    #     ##determine the distance between each edge
    #     ##once the distance is <1 exit the loop of the code 
    #         shortest_distance_to_atom = 1000 #setting arbritraury distance that much greater than the radius of all atoms 


    #         ##the plane should be the one defined by the triangle...
    #         for edge, (edge_vector, point_on_edge) in edges_and_points.items():
    #             #triangle_atom_assignment_upper[triangle_idex_num] = []
    #             triangle_atom_assignment_upper[triangle_idex_num] = 0  

    #             if shortest_distance_to_atom > atm_radius:
    #                 #vector from point of vector to coords_atom_plane_of_triangle
                    
    #                 vector_to_atom = atom_coords- point_on_edge
    #                 area_parallelogram = np.linalg.norm(np.cross(vector_to_atom, edge_vector))  #this is giving me coordinates and not a single value....
    #                 base_parallelogram = np.linalg.norm(edge_vector)

    #                 shortest_distance_to_atom = np.divide(area_parallelogram,base_parallelogram)
    #             if shortest_distance_to_atom <= atm_radius:
    #                 if atm_type == "n":  #add an if statement where if the atom is polar exit and move onto a new triangle
    #                     triangle_atom_assignment_upper[triangle_idex_num] += 1 
    #                     # if triangle_idex_num not in triangle_atom_assignment_upper:
    #                     #     triangle_atom_assignment_upper[triangle_idex_num] = []
    #                     #triangle_atom_assignment_upper[triangle_idex_num].append((atom_id).strip()) 
    #                     #triangle_atom_assignment_upper[triangle_idex_num].append((atom_id).strip()) 
    #                     ##rather than appends the atom id, maybe just add 1 
    #                     stop_triangle = True
    #                     break 
    #                 else:
    #                     triangle_atom_assignment_upper[triangle_idex_num] += .001
    #                     # if triangle_idex_num not in triangle_atom_assignment_upper:
    #                     #     triangle_atom_assignment_upper[triangle_idex_num] = []
    #                     #triangle_atom_assignment_upper[triangle_idex_num].append((atom_id).strip()) 
    #                     break
             
    # # Save as a pickle file
    # with open("triangle_atom_assignment_upper_test3.pickle", "wb") as f:
    #     pickle.dump(triangle_atom_assignment_upper, f) 




