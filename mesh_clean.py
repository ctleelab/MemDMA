import pickle 
from pathlib import Path
import numpy as np
import MDAnalysis as mda

#loading data from mc and formating data for analysis
def format_data(X,Y,Zupper, Zlower):
    #remove the extra dimension from Z_upper
    #(X,Y) should be vertices of the quadrilateral

    Zupper_2d = np.divide(Zupper.squeeze(), 10) #removes the dimension representing the Z column
    Zlower_2d = np.divide(Zlower.squeeze(), 10) #removes the dimension representing the Z column
    #Y_np = np.divide(np.array(Y), 10) #rescaled to be in nm 
    Y_np = np.array(Y) 
    X_np = np.array(X)
    #np.savetxt("Y_np", Y_np)
    #X_np = np.divide(np.array(X),10) #rescaled to be in nm 
    #np.savetxt("X_np", X_np)

    return X_np, Y_np, Zupper_2d, Zlower_2d

#determines which residues belong to the upper vs lower leaflet
def residue_leaflet_determine(u, selection):
    import leaflet_finder
    import MDAnalysis as mda 
    leaflet, residue_leaflet = leaflet_finder.determine_leaflets(u, selection)

    return leaflet, residue_leaflet

#leaflet, residue_leaflet = residue_leaflet_determine(u,selection)

#{triangle number : centroid [(x1,y1,z1), normal vector (x2,y2,z2)}
#goal: calculate the center for the each triangle
def normal_vector_triangles(X,Y,Zupper_2d, Zlower_2d):

    #center of triangle = (Ax + Bx + Cx)/3, (Ay + By + Cy)/3, (Az + Bz + Cz)/3
    centroid_triangles_upper = {}
    centroid_triangles_lower = {}
    normal_vector_triangles_upper = {}
    normal_vector_triangles_lower = {}
    coords_triangles_upper = {}
    coords_triangles_lower = {}
    edges_triangles_upper = {}
    edges_triangles_lower = {}
    quadrilateral_upper = {}
    quadrilateral_lower = {} 
    # vertices_quadrilateral = set()
    # connections_quadrilateral = []
    nx, ny = X.shape
    triangle_upper = 1
    triangle_lower = 1 
    quadrilateral_number = 1
    print("This is nx:", nx)
    print("This is ny:", ny)
    for i in range(nx-1):
        for j in range(ny-1):
            #calculate the centroids of the triangles by breaking down each quadrilateral into two triangles
            #index of the points in the quadrilateral
            p0 = (i,j)
            p1 = (i+1,j)
            p2 = (i+1, j+1)
            p3 = (i, j+1)

            #values of the (X,Y,Z) for each point that composes a quadrilateral
            quad_upper = np.array([[X[p0], Y[p0], Zupper_2d[p0]],
                                   [X[p1], Y[p1], Zupper_2d[p1]],
                                   [X[p2], Y[p2], Zupper_2d[p2]],
                                   [X[p3], Y[p3], Zupper_2d[p3]]]
                                   )
            # for point in quad_upper:
            #     vertices_quadrilateral.add(tuple(point))
            # connections_quadrilateral.append([tuple(quad_upper[0]), tuple(quad_upper[1]), tuple(quad_upper[2]), tuple(quad_upper[3])])
            
            #triangle 1: [p0,p1,p2]
            traingle1_xcord = (quad_upper[0,0] + quad_upper[1,0] + quad_upper[2,0])/3 
            traingle1_ycord = (quad_upper[0,1] + quad_upper[1,1] + quad_upper[2,1])/3 
            traingle1_zcord = (quad_upper[0,2] + quad_upper[1,2] + quad_upper[2,2])/3 
            triangle1_centroid = (traingle1_xcord, traingle1_ycord, traingle1_zcord)
            centroid_triangles_upper[triangle_upper] = np.array(triangle1_centroid)

            edge1_traingle1 = quad_upper[1] - quad_upper[0]
            edge2_traingle1 = quad_upper[2] - quad_upper[1]
            edge3_traingle1 = quad_upper[0] - quad_upper[2]


            coords_triangles_upper[triangle_upper] = [[X[p0], Y[p0], Zupper_2d[p0]],[X[p1], Y[p1], Zupper_2d[p1]], [X[p2], Y[p2], Zupper_2d[p2]]]

            edges_triangles_upper[triangle_upper] = [edge1_traingle1, edge2_traingle1, edge3_traingle1]

            #calculating normal vectors to each triangle
            #upper, triangle 1: will be vector from p0 to p1 and p1 to p2
            v1_1 = quad_upper[1] - quad_upper[0]
            v1_2 = quad_upper[2] - quad_upper[0]

            #cross product between the two vector
            normal1 = np.cross(v1_1, v1_2)
            #normal1_unit =  np.linalg.norm(normal1) 
            normal_vector_triangles_upper[triangle_upper]= normal1
            
            # #assigning the quadrilateral number to each triangle
            # if quadrilateral_number not in quadrilateral_upper:
            #     quadrilateral_upper[quadrilateral_number] = []
            # quadrilateral_upper[quadrilateral_number].append(triangle_upper)

            #increase key number of triangle
            triangle_upper += 1 


            #triangle 2: [p0,p2,p3]
            traingle2_xcord = (quad_upper[0,0] + quad_upper[2,0] +  quad_upper[3,0])/3 
            traingle2_ycord = (quad_upper[0,1] + quad_upper[2,1] + quad_upper[3,1])/3 
            traingle2_zcord = (quad_upper[0,2] + quad_upper[2,2] + quad_upper[3,2])/3 
            triangle2_centroid = (traingle2_xcord, traingle2_ycord, traingle2_zcord)
            centroid_triangles_upper[triangle_upper] = np.array(triangle2_centroid)

            edge1_traingle2 = quad_upper[2] - quad_upper[0]
            edge2_traingle2 = quad_upper[3] - quad_upper[2]
            edge3_traingle2 = quad_upper[0] - quad_upper[3]

            coords_triangles_upper[triangle_upper] = [[X[p0], Y[p0], Zupper_2d[p0]],[X[p2], Y[p2], Zupper_2d[p2]], [X[p3], Y[p3], Zupper_2d[p3]]] 

            edges_triangles_upper[triangle_upper] = [edge1_traingle2, edge2_traingle2, edge3_traingle2]

            #upper, triangle 2: will be vector from p2 to p3 and p3 to p1
            v2_1 = quad_upper[3] - quad_upper[2]
            v2_2 = quad_upper[0] - quad_upper[3]
            normal2 = np.cross(v2_1, v2_2)

            #normal2_unit =  np.linalg.norm(normal2) 
            normal_vector_triangles_upper[triangle_upper]= normal2

            #averagenormal
            average_normal_upper = np.average([normal1, normal2], axis = 0)

            #assigning the quadrilateral number to each triangle
            if quadrilateral_number not in quadrilateral_upper:
                quadrilateral_upper[quadrilateral_number] = []
            quadrilateral_upper[quadrilateral_number].append(average_normal_upper)      

            triangle_upper += 1

            ###############################
            #####start of lower leaflet####
            ###############################

            quad_lower = np.array([[X[p0], Y[p0], Zlower_2d[p0]],
                                   [X[p1], Y[p1], Zlower_2d[p1]],
                                   [X[p2], Y[p2], Zlower_2d[p2]],
                                   [X[p3], Y[p3], Zlower_2d[p3]]]
                                   )
            

            #triangle 3: [p0,p1,p2]
            traingle3_xcord = (quad_lower[0,0] + quad_lower[1,0] + quad_lower[2,0])/3 
            traingle3_ycord = (quad_lower[0,1] + quad_lower[1,1] + quad_lower[2,1])/3 
            traingle3_zcord = (quad_lower[0,2] + quad_lower[1,2] + quad_lower[2,2])/3 
            triangle3_centroid = (traingle3_xcord, traingle3_ycord, traingle3_zcord)
            centroid_triangles_lower[triangle_lower] = np.array(triangle3_centroid)


            edge1_traingle3 = quad_lower[1] - quad_lower[0]
            edge2_traingle3 = quad_lower[2] - quad_lower[1]
            edge3_traingle3 = quad_lower[0] - quad_lower[2]

            coords_triangles_lower[triangle_lower] = [[X[p0], Y[p0], Zlower_2d[p0]],[X[p1], Y[p1], Zlower_2d[p1]], [X[p2], Y[p2], Zlower_2d[p2]]]


            edges_triangles_lower[triangle_lower] = [edge1_traingle3, edge2_traingle3, edge3_traingle3]

            #lower, triangle 3: will be vector from p0 to p1 and p1 to p2
            v3_1 = quad_lower[1] - quad_lower[0]
            v3_2 = quad_lower[2] - quad_lower[0]
            normal3 = np.cross(v3_1, v3_2)
            #normal3_unit =  np.linalg.norm(normal3) 
            normal_vector_triangles_lower[triangle_lower]= normal3


            # #assigning the quadrilateral number to each triangle
            # if quadrilateral_number not in quadrilateral_lower:
            #     quadrilateral_lower[quadrilateral_number] = []
            # quadrilateral_lower[quadrilateral_number].append(triangle_lower)

            triangle_lower += 1

            #triangle 4: [p0,p2,p3]
            traingle4_xcord = (quad_lower[0,0] + quad_lower[1,0] + quad_lower[2,0])/3 
            traingle4_ycord = (quad_lower[0,1] + quad_lower[1,1] + quad_lower[2,1])/3 
            traingle4_zcord = (quad_lower[0,2] + quad_lower[1,2] + quad_lower[2,2])/3 
            triangle4_centroid = (traingle4_xcord, traingle4_ycord, traingle4_zcord)
            centroid_triangles_lower[triangle_lower] = np.array(triangle4_centroid)

            edge1_traingle4 = quad_lower[2] - quad_lower[0]
            edge2_traingle4 = quad_lower[3] - quad_lower[2]
            edge3_traingle4 = quad_lower[0] - quad_lower[3]

            coords_triangles_lower[triangle_lower] = [[X[p0], Y[p0], Zlower_2d[p0]],[X[p2], Y[p2], Zlower_2d[p2]], [X[p3], Y[p3], Zlower_2d[p3]]]

            edges_triangles_lower[triangle_lower] = [edge1_traingle4, edge2_traingle4, edge3_traingle4]

            #lower, triangle 4: will be vector from p2 to p3 and p3 to p1
            v4_1 = quad_lower[3] - quad_lower[2]
            v4_2 = quad_lower[0] - quad_lower[3]
            normal4 = np.cross(v4_1, v4_2)
            #normal4_unit =  np.linalg.norm(normal4) 
            normal_vector_triangles_lower[triangle_lower] = normal4
            
            #averagenormal between the two triangles that make up the quadrilateral
            average_normal_lower = np.average([normal3, normal4], axis = 0)
            
            #assigning the quadrilateral number to each triangle
            if quadrilateral_number not in quadrilateral_upper:
                quadrilateral_upper[quadrilateral_number] = []
            quadrilateral_upper[quadrilateral_number].append(average_normal_lower)

            

            triangle_lower += 1 

            quadrilateral_number += 1


    # centroid_triangles_upper_np = np.array(centroid_triangles_upper)
    # centroid_triangles_upper_flatten = centroid_triangles_upper_np.reshape(centroid_triangles_upper_np.shape[0], -1)
    # np.savetxt("centroid_triangles_upper", centroid_triangles_upper_flatten)

    # normal_vector_triangles_upper_np = np.array(normal_vector_triangles_upper)
    # normal_vector_triangles_lower_np = np.array(normal_vector_triangles_lower)
    # centroid_triangles_lower_np = np.array(centroid_triangles_lower)


    return centroid_triangles_upper, centroid_triangles_lower, normal_vector_triangles_upper, normal_vector_triangles_lower, coords_triangles_upper, coords_triangles_lower, edges_triangles_upper, edges_triangles_lower, quadrilateral_upper, quadrilateral_lower

#centroid_triangles_upper_np, centroid_triangles_lower_np, normal_vector_triangles_upper_np, normal_vector_triangles_lower_np, normal_vector_triangles_upper_triangle1_np, normal_vector_triangles_upper_triangle2_np, vertices_quadrilateral_np, connections_quadrilateral_np =  normal_vector_triangles(X,Y,Zupper_2d, Zlower_2d, triangle_bins)


# #input triangle normal vector, and the quadrilateral the triangle belongs to 
# def average_quadrilateral_normal(normal_vector, quadrilateral_number):
#     #loop over the quadrilaterals 
#     #if the triangle name in it then append 






    




#determine which grid gets assigned to which lipid 
#residue_leaflet_assignment will tell if the lipid belongs to the uppper or lower leaflet
#uses euclidean distance between gl2 coordinates and centroid of the quadrilateral triangle to determine which grid is the closet to the lipid 
#output: grid mapping = gl2 coords with the coordinates of the closest triangle grid to it (should reset with each new residue)
#output: closest_mapping = residue:triangular grid with the shortest ED
def closest_grid(gl2_coords, empty_res, residue_leaflet_assignment, centroid_triangles_upper, centroid_triangles_lower):

    #create a copy of the empty_res this is where to closet grid will be stored 
    grid_mapping = empty_res.copy()
    closest_mapping = empty_res.copy()

    
    #loop through the gl2_coords
    #if the atom is a GL2 then do the following
    for residue in gl2_coords:
        residue_gl2_coords = gl2_coords[residue]
        gl2_ED = {} ##triangle number: all EDs

        if residue in residue_leaflet_assignment["upper"]:
            for triangle, centroid_coords in centroid_triangles_upper.items():
                centroid_coords = centroid_triangles_upper[triangle]
                difference = np.array(centroid_coords) - np.array(residue_gl2_coords)  ##must convert to numpy array to do the subtraction
                distance =  np.sqrt(np.sum(difference**2))
                gl2_ED[triangle] = distance
            grid_mapping[residue] = gl2_ED
            min_triangle = min(gl2_ED, key = gl2_ED.get)
            closest_mapping[residue] = min_triangle
        if residue in residue_leaflet_assignment["lower"]:
            for triangle, centroid_coords in centroid_triangles_lower.items():
                centroid_coords = centroid_triangles_upper[triangle]
                difference = np.array(centroid_coords) - np.array(residue_gl2_coords)
                distance = np.sqrt(np.sum(difference**2))
                gl2_ED[triangle] = distance
            grid_mapping[residue] = gl2_ED
            min_triangle = min(gl2_ED, key = gl2_ED.get)
            closest_mapping[residue] = min_triangle
        
        
        ##determining the lowest ED and the triangle associated with it
        ##gets the triangle with the smallest ED to the gl2 bead

        ##pair the min_triangle to it's gl2 bead
        

        # if residue in residue_leaflet_assignment["lower"]:
        #     difference = centroid_triangles_lower - coords
        #     distances = np.sqrt(np.sum(difference**2, axis = 1))
        #     ED_lower = {i: d for i, d in enumerate(distances)}
        #     min_key = min(ED_lower, key = ED_lower.get)
        #     grid_mapping[residue] = min_key

    return closest_mapping