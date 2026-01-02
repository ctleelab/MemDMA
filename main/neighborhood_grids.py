import numpy as np
from scipy.spatial import cKDTree


def cKDtree_grid_neighborhood(grid_centroids: dict):
    """
    This function creates a CKDtree from the grid_centroids.

    Args:
        grid_centroids (dict): {grid number: centroid (x,y,z)}

    Returns:
        tree: The cKDTree object for spatial queries.
        grid_keys: List of grid numbers.
        centroids: Numpy array of centroid coordinates for each grid
    """

    grid_keys = np.arange(0, len(grid_centroids))  # creates a list of the grid numbers
    centroids = np.array(
        [grid_centroids[k][:3] for k in grid_keys]
    )  # take the coordinates for each grid in grid_centroids and creates a numpy array

    tree = cKDTree(centroids)

    return tree, grid_keys, centroids


def neighborhood(cutoff: int, grid_number: int, tree, grid_keys: list, centroids):
    """
    This function creates a neighborhood for each grid with a specified cutoff value.

    Args:
        cuttoff (int): The neighborhood cutoff distance in nm.
        grid_number (int): The grid number for which to find neighborhoods.
        tree: The cKDTree object for spatial queries.
        grid_keys: List of grid numbers
        centroids: Numpy array of centroid coordinates for each grid

    Returns:
        list: A list of nearby grids determined through cDKtree query.
    """
    # create neighborhood grid for each grid, using the centroid coordinates of the grid within the distance of the cutoff value
    nearby_grids = tree.query_ball_point(centroids[grid_number], cutoff)

    # convert the nearby_grid index to a actual grid value
    nearby_keys = [grid_keys[j] for j in nearby_grids]
    return nearby_keys
