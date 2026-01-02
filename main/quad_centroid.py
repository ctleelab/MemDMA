import numpy as np
import neighborhood_grids as ng


def centroid_quadrilaterals(X, Y, Zupper_2d, Zlower_2d):
    # #trouble-shooting dimension of array
    # for key, value in mc.items():
    #     arr = np.asarray(value)
    #     print(f"{key}: shape = {arr.shape}, ndim = {arr.ndim}")
    # print(mc["Z_lower"])

    # grid:(x,y,z) for upper and lower
    # this gives each vertex of each quadrilateral composing the grid
    # ravel flattens the numpy array to 1D (ex: [1,2,3],[4,5,6] --> [1,2,3,4,5,6])
    coords_upper = np.column_stack((X.ravel(), Y.ravel(), Zupper_2d.ravel()))
    coords_lower = np.column_stack((X.ravel(), Y.ravel(), Zlower_2d.ravel()))

    # computing the center of each individual quadrilateral
    # tells the shape of the mesh and the quadrilaterals that compose it
    nrows, ncols = X.shape

    # creating empty list
    centroids_upper = []
    centroids_lower = []

    centroids_upper_dict = {}
    centroids_lower_dict = {}

    key = 0
    for i in range(nrows - 1):
        for j in range(ncols - 1):
            # calculates the centroids for each grid
            quad_upper = np.array(
                [
                    [X[i, j], Y[i, j], Zupper_2d[i, j]],
                    [X[i + 1, j], Y[i + 1, j], Zupper_2d[i + 1, j]],
                    [X[i, j + 1], Y[i, j + 1], Zupper_2d[i, j + 1]],
                    [X[i + 1, j + 1], Y[i + 1, j + 1], Zupper_2d[i + 1, j + 1]],
                ]
            )
            quad_lower = np.array(
                [
                    [X[i, j], Y[i, j], Zlower_2d[i, j]],
                    [X[i + 1, j], Y[i + 1, j], Zlower_2d[i + 1, j]],
                    [X[i, j + 1], Y[i, j + 1], Zlower_2d[i, j + 1]],
                    [X[i + 1, j + 1], Y[i + 1, j + 1], Zlower_2d[i + 1, j + 1]],
                ]
            )
            # print(f"QUAD UPPER:", quad_upper)
            # by adding axis = 0, it averaging down rows, so averaging across x, y,z independently
            mean_quad_upper = quad_upper.mean(axis=0)
            mean_quad_lower = quad_lower.mean(axis=0)
            # print(f"MEAN UPPER:", mean_quad_upper)

            centroids_upper.append(mean_quad_upper)
            centroids_lower.append(mean_quad_lower)

            # create a dictionary that stores grid: coords
            centroids_upper_dict[key] = mean_quad_upper
            centroids_lower_dict[key] = mean_quad_lower

            key += 1

    # covert back to numpy array
    centroids_upper_matrix = np.array(centroids_upper)
    centroids_lower_matrix = np.array(centroids_lower)

    return (
        centroids_upper_matrix,
        centroids_lower_matrix,
        centroids_upper_dict,
        centroids_lower_dict,
    )


# get the index of the smallest value (by using the function np.argmin)
# the output of this is the same for all keys
# grid_dict: grid with it's neighborhood
# i
def grid_assignment(
    gl2_coords: dict,
    grid_matrix,
    tree_upper,
    grid_keys_upper: list,
    centroids_upper,
    cutoff: int,
    atom_features: dict,
):
    """
    Assigns a grid to a GL2 bead then computes its neighborhood using cKDTree.
    Args:
        gl2_coords (dict): A dictionary containing GL2 bead coordinates.
        grid_matrix (np.ndarray): A matrix of grid coordinates.
        tree_upper (cKDTree): A cKDTree object for efficient spatial queries.
        grid_keys_upper (list): A list of grid keys corresponding to the grid matrix.
        centroids_upper (np.ndarray): An array of centroid coordinates for the upper grid.
        cutoff (int): The cutoff distance for neighborhood determination.

    Returns:
        gl2_grid_assignment (dict): A dictionary with GL2 (key): (closet grid_index, grid_coords of closest grid_index, neighborhood_list of the grid_idex).
        lipid_grid_assignment (dict): A dictionary with lipid (key): (neighborhood_list of the grid_idex)
    """

    gl2_grid_assignment = {}  # create an empty dictionary with the same keys
    lipid_grid_assignment = {}

    for key, gl2_points in gl2_coords.items():
        # print(grid_matrix.shape)
        distance = np.linalg.norm(
            (grid_matrix - gl2_points), axis=1
        )  # compute the euclidean distance between grid and gl2, and sum along each row
        min_distance_index = np.argmin(
            distance
        )  # return the smallest value of the distance array index value

        grid_coords = grid_matrix[min_distance_index]
        # rather then grid_coords I want grid_number

        # match the grid_coords to the grid_dict and pull out the key
        list_neighborhood = ng.neighborhood(
            cutoff, min_distance_index, tree_upper, grid_keys_upper, centroids_upper
        )

        gl2_grid_assignment[key] = (min_distance_index, grid_coords, list_neighborhood)

        lipid_res = atom_features[key][4]

        # store lipid_res as the key and the neighborhood list as the value
        lipid_grid_assignment[lipid_res] = list_neighborhood

    for k, v in list(gl2_coords.items())[:5]:
        print(k, np.round(v, 3))

    return gl2_grid_assignment, lipid_grid_assignment


# gl2_neighborhood_grids: {gl2_atom_number: grid_number, (coords of grid), list of neighborhood grids}
# all_atom: {atom_number: resnumber, coords}
# gl2_beads_resnumber: {gl2_atom_number: resnumber}


def atom_grid_mapping(
    all_atom: dict, gl2_neighborhood_grids: dict, gl2_beads_resnumber: dict
):
    """
    Append a list of neighborhood grids to the all_atoms dictionary.
    Args:
        all_atom (dict): {atom number: res_number, coords}
        gl2_neighborhood_grids (dict): {A dictionary with GL2 atom number: (closet grid_index, grid_coords of closest grid_index, neighborhood_list of the grid_idex)}
        gl2_beads_resnumber (dict): {atom number of gl2 bead: res_number}
    Returns:
        all_atom (dict): {atom number: res_number, coords, list of neighborhood grid number}
    """

    # Precompute reverse lookup. Mapping residue number -> gl2_atom
    res_to_gl2_atom = {v: str(k).strip() for k, v in gl2_beads_resnumber.items()}

    for key1, v1 in all_atom.items():

        res_number = v1[0]  # pulls out residue number

        # Find corresponding GL2 atom for each residue number
        gl2_atom = res_to_gl2_atom.get(res_number)

        if (
            gl2_atom is None
        ):  # Error handeling of residues that are not present in GL2 mapping
            print("Error, residue not found in the Gl2 mapping.")
            break

        values_1 = gl2_neighborhood_grids.get(
            gl2_atom
        )  # pulls out the values of the gl2_neighborhood_grids dictionary

        neighborhood_list = values_1[
            2
        ]  # pulls out the neighborhood grid list of the gl2_neighborhood_grids dictionary

        # Append neighborhood list to atom entry and preserves the original data
        all_atom[key1] = (*v1, neighborhood_list)

    return all_atom
