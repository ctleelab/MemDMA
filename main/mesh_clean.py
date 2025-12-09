import pickle 
from pathlib import Path
import numpy as np
import MDAnalysis as mda 
import csv

#loading data from mc and formating data for analysis
def format_data(X,Y,Zupper, Zlower):
    #remove the extra dimension from Z_upper
    #(X,Y) should be vertices of the quadrilateral

    # Zupper_2d = np.divide(Zupper.squeeze(), 10) #removes the dimension representing the Z column
    # Zlower_2d = np.divide(Zlower.squeeze(), 10) #removes the dimension representing the Z column
    Zupper_2d = Zupper.squeeze() #removes the dimension representing the Z column
    Zlower_2d = Zlower.squeeze() #removes the dimension representing the Z column
    Y_np = np.array(Y) 
    X_np = np.array(X)
    #np.savetxt("Y_np", Y_np)
    #np.savetxt("X_np", X_np)

    return X_np, Y_np, Zupper_2d, Zlower_2d

#determines which residues belong to the upper vs lower leaflet
def residue_leaflet_determine(u, selection):
    import main.leaflet_finder as leaflet_finder
    import MDAnalysis as mda 
    leaflet, residue_leaflet = leaflet_finder.determine_leaflets(u, selection)

    return leaflet, residue_leaflet

#leaflet, residue_leaflet = residue_leaflet_determine(u,selection)

#{triangle number : centroid [(x1,y1,z1), normal vector (x2,y2,z2)}
#goal: calculate the center for the each triangle

def normal_vector_triangles(X, Y, Zupper_2d, Zlower_2d, normalize=True):
    nx, ny = X.shape
    n_quads = (nx - 1) * (ny - 1)
    n_tri = 2 * n_quads  # two triangles per quad per leaflet

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
    edges_triangles_upper_2vertices = {}
    edges_triangles_lower_2vertices = {}

    centroid_upper_list = []
    centroid_lower_list = []

    triangle_upper = 1
    triangle_lower = 1
    quad_num = 1

    for i in range(nx - 1):
        for j in range(ny - 1):
            # corners of the original grid
            p0 = (i, j)
            p1 = (i + 1, j)
            p2 = (i + 1, j + 1)
            p3 = (i, j + 1)

            # ---- upper leaflet ----
            quad_u = np.array([
                [X[p0], Y[p0], Zupper_2d[p0]],
                [X[p1], Y[p1], Zupper_2d[p1]],
                [X[p2], Y[p2], Zupper_2d[p2]],
                [X[p3], Y[p3], Zupper_2d[p3]],
            ])

            # two triangles: [p0,p1,p2] and [p0,p2,p3]
            tri1 = quad_u[[0, 1, 2]]
            tri2 = quad_u[[0, 2, 3]]

            # centroids
            cent1 = tri1.mean(axis=0) #average across coordinates of the triangle
            cent2 = tri2.mean(axis=0)

            centroid_triangles_upper[triangle_upper] = cent1
            centroid_triangles_upper[triangle_upper + 1] = cent2
            centroid_upper_list.extend([cent1, cent2])

            # normals
            n1 = np.cross(tri1[1] - tri1[0], tri1[2] - tri1[0]) #cross product between two edge vectors 
            n2 = np.cross(tri2[1] - tri2[0], tri2[2] - tri2[0])
            if normalize:
                n1 /= np.linalg.norm(n1) 
                n2 /= np.linalg.norm(n2) 

            normal_vector_triangles_upper[triangle_upper] = n1
            normal_vector_triangles_upper[triangle_upper + 1] = n2

            coords_triangles_upper[triangle_upper] = tri1
            coords_triangles_upper[triangle_upper + 1] = tri2

            edges_triangles_upper[triangle_upper] = [tri1[1] - tri1[0], tri1[2] - tri1[1], tri1[0] - tri1[2]]
            edges_triangles_upper[triangle_upper + 1] = [tri2[1] - tri2[0], tri2[2] - tri2[1], tri2[0] - tri2[2]]

            edges_triangles_upper_2vertices[triangle_upper] = [[tri1[0],tri1[1]], [tri1[1], tri1[2]], [tri1[0],tri1[2]]]
            edges_triangles_upper_2vertices[triangle_upper + 1] = [[tri2[0],tri2[1]], [tri2[1], tri2[2]], [tri2[0],tri2[2]]]

            avg_normal_upper = (n1 + n2) / 2
            quadrilateral_upper[quad_num] = [avg_normal_upper]

            triangle_upper += 2

            # ---- lower leaflet ----
            quad_l = np.array([
                [X[p0], Y[p0], Zlower_2d[p0]],
                [X[p1], Y[p1], Zlower_2d[p1]],
                [X[p2], Y[p2], Zlower_2d[p2]],
                [X[p3], Y[p3], Zlower_2d[p3]],
            ])

            tri3 = quad_l[[0, 1, 2]]
            tri4 = quad_l[[0, 2, 3]]
            cent3 = tri3.mean(axis=0)
            cent4 = tri4.mean(axis=0)
            centroid_triangles_lower[triangle_lower] = cent3
            centroid_triangles_lower[triangle_lower + 1] = cent4
            centroid_lower_list.extend([cent3, cent4])

            n3 = np.cross(tri3[1] - tri3[0], tri3[2] - tri3[0])
            n4 = np.cross(tri4[1] - tri4[0], tri4[2] - tri4[0])
            if normalize:
                n3 /= np.linalg.norm(n3) 
                n4 /= np.linalg.norm(n4)

            normal_vector_triangles_lower[triangle_lower] = n3
            normal_vector_triangles_lower[triangle_lower + 1] = n4

            coords_triangles_lower[triangle_lower] = tri3
            coords_triangles_lower[triangle_lower + 1] = tri4

            edges_triangles_lower[triangle_lower] = [tri3[1] - tri3[0], tri3[2] - tri3[1], tri3[0] - tri3[2]]
            edges_triangles_lower[triangle_lower + 1] = [tri4[1] - tri4[0], tri4[2] - tri4[1], tri4[0] - tri4[2]]

            edges_triangles_lower_2vertices[triangle_upper] = [[tri3[0],tri3[1]], [tri3[1], tri3[2]], [tri3[0],tri3[2]]]
            edges_triangles_lower_2vertices[triangle_upper + 1] = [[tri4[0],tri4[1]], [tri4[1], tri4[2]], [tri4[0],tri4[2]]]

            avg_normal_lower = (n3 + n4) / 2
            quadrilateral_lower[quad_num] = [avg_normal_lower]

            triangle_lower += 2
            quad_num += 1
    #####################
    ###troubleshooting###
    #####################
    
    # z_coords = []
    # for arr in centroid_upper_list:
    #     z_coord = (arr[2])
    #     z_coords.append(z_coord)
    # print("min:", min(z_coords))
    # print("max:", max(z_coords))

    centroid_triangles_upper_matrix = np.vstack(centroid_upper_list)
    centroid_triangles_lower_matrix = np.vstack(centroid_lower_list)


    with open(
        f"centroid_triangle.csv",
        "w",
        newline="",
    ) as f:
        writer = csv.writer(f)
        writer.writerow(["triangle_id", "coordinates"])  # optional header
        for key, value in centroid_triangles_upper.items():
            writer.writerow([key, value])
    return (
        centroid_triangles_upper,
        centroid_triangles_lower,
        normal_vector_triangles_upper,
        normal_vector_triangles_lower,
        coords_triangles_upper,
        coords_triangles_lower,
        edges_triangles_upper,
        edges_triangles_lower,
        quadrilateral_upper,
        quadrilateral_lower,
        centroid_triangles_upper_matrix,
        centroid_triangles_lower_matrix,
        edges_triangles_upper_2vertices,
        edges_triangles_lower_2vertices
    )


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