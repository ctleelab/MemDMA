#goal: this will be the main script in which all functions will be runned 
#notes:
###:denotes plan of code
##: this denotes comments by Cailyn that are made to understand the fundamentals of general python programming
#: this denotes comments essential for understanding the code


import argparse
import sys
import pdb as pdb
import numpy as np
import re 
import MDAnalysis as mda 

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

##this chunk of code will only be run if the script is executed directly
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
        ##will output 3 dictionaries containing 1) taking lipid_bead head name: lipid name(3-letters) 2) lipid name: GL2 bead head  3)lipid name (3-letters)
        DICT_3L, RESNAME_GLYC, LIPID = p.set_params(args.paramFile)

        print("RESNAME_GLYC:", RESNAME_GLYC)


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


    ##tells you the column at which the residue number is stored for .gro files
    startID = 0
    endID = 5
    

    ###############################################################
    #####Building a dictionary {lipidname: [residue numbers]} #####
    ###############################################################
    zaxis_byatoms = []
    xaxis_byatoms = []
    yaxis_byatoms = []
    dicoMb = {}
    ions = ["W", "ION"]
    residue_dict = {} 
    residue_coords_gl2 = {}


    for i, line_atom in enumerate(grolines):
        if i > 1 and  i < (len(grolines)-1): ##skips the first two lines of the .gro file

            #indexing like this in a gro file is not consistent, as values keeps increasing it changes the length of the index
            #will seperate based on the numbers_text
            res_name_number=line_atom[0:12].strip() ##grabs residue name with lipid name and takes out white sapce
            match = re.match(r"(\d+)([A-Za-z]+)", res_name_number) ##\d matches any digit, + means one or more digits/letters, matches letters A-Za-z
            if match:
                res_number = int(match.group(1)) ##grabs residue number
                res_name = str(match.group(2))

            #list of atom types 
            list_types_atoms = util.list_types_atoms

            atm_name_number=str(line_atom[12:22].strip()) ##grabs the atom name and atom number and takes out whitespace
            
            ##splits the atom name from atom number 
            for j in list_types_atoms:
                if atm_name_number.startswith(j): 
                    atm_name = j
                    atom_number = atm_name_number[len(j):]
                    break 

            x_coord = float(line_atom[20:28].strip())##grabs the x coordinate
            y_coord = float(line_atom[28:36].strip()) ##grabs the y coordinate
            z_coord = float(line_atom[36:44].strip()) ##grabs the z coordinate
            xaxis_byatoms.append(x_coord)  
            yaxis_byatoms.append(y_coord) 
            zaxis_byatoms.append(z_coord)  

            key = res_name+  "_" + atm_name ##combines resname with atom name (ex: DOP_NH3)

            #atm_key = res_name + " " + atm_

            # ##I don't know if I need this code to create the DICT_3L...since I am working with .gro files
            # if res_name in LIPID: 
            #     if key in DICT_3L:
            #         if DICT_3L[key] in dicoMb: ##finds the three letter code for atom, only does it for the head group, then stores the residue number
            #             dicoMb[DICT_3L[key]].append(res_number) ##stores the residue number to the resname_atomname
            #             #print(key)
            #             #print(int(line_atom[startID:endID].strip()))                    
            #         else:
            #             dicoMb[DICT_3L[key]] = []
            #             dicoMb[DICT_3L[key]].append(res_number)

            ##creates a residue_dict {residue number: []}
            ##creates a residue_dict {residue number: [x,y,z]}
            if res_name in LIPID:
                if res_name in RESNAME_GLYC:
                    if atm_name == RESNAME_GLYC[res_name]: 
                        if res_number not in residue_dict:
                            if res_number not in residue_coords_gl2:
                                residue_dict[res_number] = []
                                residue_coords_gl2[res_number] = []
                                residue_coords_gl2[res_number].append(x_coord)
                                residue_coords_gl2[res_number].append(y_coord)
                                residue_coords_gl2[res_number].append(z_coord)
                    
            if res_name in ions:
                continue

    
    ##################################################################
    ######CHUNK OF CODE THAT IS NOT NEEDED SINCE USING GRO FILES######
    ##################################################################

    # ###Checking to ensure that the total residue number does not exceed 9999 (PDB format limit) 
    # nb_res=0
    # ##goes through each key in the dicionary and adds it count to the nb_res
    # for key in dicoMb:
    #     nb_res+=len(dicoMb[key])
    # if nb_res > 9999:
    #     print("The number of lipids in your membrane overtakes the PDB format (>9999)!!")
    #     print("This point is not supported by PackMem, this tool may not work properly !")       
    # ##Checks to see if there are duplicate in the residue 
    # flD=d.detect_duplicate(dicoMb)
    # if flD:
    #     print("You have duplicate in lipids residue_number")    
    #     print("This point is not supported by PackMem, this tool may not work properly !") 
    #     print("Please, renumerate your residue_number to avoid this problem (for example using editconf Gromacs tool)")  
    #     sys.exit()
    
    # ##Determine min and max values of the coordinates (x,y,z) that we extracted
    xmin, xmax, xmean = l.min_max(xaxis_byatoms)
    ymin, ymax, ymean = l.min_max(yaxis_byatoms)
    zmin, zmax, zmean = l.min_max(zaxis_byatoms)
    
    ##################################################################
    ######CHUNK OF CODE THAT IS NOT NEEDED SINCE USING GRO FILES######
    ##################################################################
    # ##modify the lipid name in data (3L->new3L to take account problem POPvsPOPS/POPE)
    # ##DOP --> DOE 
    # ##can tell the difference between DOPC vs DOPE lipids by DOPC(NC3) DOPE(NH3)s
    
    # ##checking the format of grolines before and after modification 
    # ##Note: in the future might want to try run PackMem with .gro files
    # print("old pbd file:", grolines[10:20])
    # grolines = pdb.modifyPDBdata(grolines, dicoMb, startID, endID)
    # print("new pdb file:", grolines[10:20])


    ################################################
    ########## Start of the real code ##############
    ################################################

    ###Here PackMem would normaly determine the upper and lower leaflets and slices based 
    ###take GL2 bead determines the grid each bead in the lipid will fall into, and go from there 
    
    ###quadrilateral mesh to create the quadrilateral mesh 
    #right now it takes in .gro, selection, path variables
    #I will be executing the quadrilateral mesh function here to get the mesh grid

    #inputs: grofile, selection of beads for mesh, step size for mesh grid (nm)
    #outputs: (X, Y, Z_upper, Z_lower) 
    ##note this also calcualtes the mesh information
    X, Y, zupper, zlower, leaflet, residue_leaflet_old, leaflet_residues_new, z_midplane = qm.mesh(args.filename, "name GL2", 1)

    ##formating data for analysis
    X_f, Y_f, zupper_f, zlower_f = m.format_data(X,Y, zupper,zlower)

    ###turns the quadrilateral mesh into a triangular mesh and computes the normal vector for each triangle in the mesh
    ###each key represents the triangle grid number
    ###normal_vector_triangles: {triangle number: (x,y,z) of normal vector}
    ###centroid_triangles: {triangle number: (x,y,z) of the centroid of the triangle}
    centroid_triangles_u, centroid_triangles_l, norm_vector_triangles_up, norm_vector_triangles_l, coords_triangles_upper, coords_triangles_lower = m.normal_vector_triangles(X_f,Y_f,zupper_f,zlower_f)

    ##############################
    ###grid assignment per atom###
    ##############################

    ##iterate over the triangles
    ##iterate over the atoms 
    ##project the bead onto the plane

    atom_ray_tracing_upper = {}
    atom_ray_tracing_lower = {}

    for triangle_num in centroid_triangles_u: 
        #iterates over the atoms
        print(triangle_num)
        #print(len(centroid_triangles_u)) #39156
        atom_features = pp.atom_features(grolines, list_types_atoms,leaflet_residues_new)
        for atom_number in atom_features:
            if atom_features[atom_number][1] == "upper": ##check the value of the dictionary 

                ##pull all necessary information 
                normal_vector_triangle = np.array(norm_vector_triangles_up[triangle_num]) #grabs the normal vector for the current triangle grid
                centroid_triangle = np.array(centroid_triangles_u[triangle_num])

                #project point onto plane, using the norm_vector_triangle as the direction vector
                w_vector = np.subtract(centroid_triangle, atom_features[atom_number][2]) #vector from centroid of triangle to atom coord 
                scalar = np.dot(w_vector,normal_vector_triangle)
                n_contribution_in_w_vector = np.multiply(scalar, normal_vector_triangle)
                vector_orthogonal_normal = np.subtract(w_vector, n_contribution_in_w_vector)


                #for each trinangle I will have coords for each atom
                atom_ray_tracing_upper[triangle_num] = {}
                atom_ray_tracing_upper[triangle_num][atom_number] = vector_orthogonal_normal  #atom point projected onto traingle plane

            if leaflet == "lower":
                continue


    print("done with upper leaflet")

    for triangle_num in centroid_triangles_l: 
        #iterates over the atoms
        print(triangle_num)
        #print(len(centroid_triangles_u)) #39156
        atom_features = pp.atom_features(grolines, list_types_atoms,leaflet_residues_new)
        for atom_number in atom_features:
            if atom_features[atom_number][1] == "lower": ##check the value of the dictionary 

                ##pull all necessary information 
                normal_vector_triangle = np.array(norm_vector_triangles_l[triangle_num]) #grabs the normal vector for the current triangle grid
                centroid_triangle = np.array(centroid_triangles_l[triangle_num])

                #project point onto plane, using the norm_vector_triangle as the direction vector
                w_vector = np.subtract(centroid_triangle, atom_features[atom_number][2]) #vector from centroid of triangle to atom coord 
                scalar = np.dot(w_vector,normal_vector_triangle)
                n_contribution_in_w_vector = np.multiply(scalar, normal_vector_triangle)
                vector_orthogonal_normal = np.subtract(w_vector, n_contribution_in_w_vector)


                #for each trinangle I will have coords for each atom
                atom_ray_tracing_lower[triangle_num] = {}
                atom_ray_tracing_lower[triangle_num][atom_number] = vector_orthogonal_normal  #atom point projected onto traingle plane

            if leaflet == "upper":
                continue


    #################################################
    ##now determine what atoms are in each triangle##
    #################################################
    ###outline of code###
    ##1) calculate the shortest distance from from point to edge of the triangle
        ##will have to do for each edge 
    
    # triangle_atom_assignment = {}
    # for triangle in atom_ray_tracing_upper:
    #     for each triangle loop through atom_number:
    #         for each edge in triangle: #determine the vectors for the edges in the mesh_clean function? 
    #             edge1 = P0 - P1
    #             edge2 = P1 - P2
    #             edge3=  P2 - P0 
    #             shortest_distance = np.divide(np.absolute(np.cross(vector1, vector2)),np.absolute(vector1)): 
    #             if shortest_distance < radius_of_atom:
    #                 then add it to the triangle list:
    #                     for triangle_atom_assignment[] = atom_number

    
        
        


    # #####################
    # ###grid assignment###    
    # #####################
    # #next determine for each lipid which grid is the closest, then use that to iterate up and down the z-axis 
    # #output: dictionary{lipid number: index of the grid number it's associated with}


    # #input: GL2 bead and its coordinates, centroids of triangle, empty list of residue numbers, upper vs lower leaflet assignment based on residue number
    # #then calculate ED to determine which grid is the closest for that lipid
    # #output: closest_mapping = {residue number:triangle number} 
    # closest_mapping = m.closest_grid(residue_coords_gl2,residue_dict,leaflet_residues_new, centroid_triangles_u, centroid_triangles_l)
    # print("completed closest_mapping")
  
    
    # # #input: gro file, residue assignment to leaflet, grid_mapping, vectors 
    # # #output: atm_grid {atom number: (x,y,z) coordinates of the normal vector of the grid that is closest to the lipid that the atom belongs to}
    # atm_normal_vectors = gs.atm_normal_vector_assignment(grolines, leaflet_residues_new, closest_mapping, norm_vector_triangles_up, norm_vector_triangles_l)
    # print("completed atm_normal_vectors")


    # atm_centroids = gs.atm_centroid_assignment(grolines, leaflet_residues_new, closest_mapping, centroid_triangles_u, centroid_triangles_l)
    # print("completed atm_centroids")

    # #solve for intersection of z_midplane and normal vectors
    # atm_lower_limit = gs.atm_lower_limit_assignment(grolines,leaflet_residues_new, atm_normal_vectors, atm_centroids)


    # #output: {atom_number: [centroid of all slices (x,y,z)]}
    # upper_slices = gs.grid_slices(grolines, leaflet_residues_new, atm_centroids, atm_normal_vectors, f_mid)
    
    # first_key = list(upper_slices.keys())[0]
    # print("Key:", first_key)
    # print("value of key:", upper_slices[first_key])
    






