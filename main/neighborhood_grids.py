import numpy as np
from scipy.spatial import cKDTree


def cKDtree_grid_neighborhood(grid_centroids:dict):
    """
    This function creates a CKDtree from the grid_centroids.

    Args:
        grid_centroids (dict): A dictionary with grid numbers as keys and their centroids as values.

    Returns:
        tree: The cKDTree object for spatial queries.
        grid_keys: List of grid numbers.
        centroids: Numpy array of centroid coordinates for each grid 
    """

    grid_keys = list(grid_centroids.keys())  #creates a list of the keys of all centroids of the grid points
    centroids = np.array([grid_centroids[k][:3] for k in grid_keys]) #take the coordinates for each grid in grid_keys and creates a numpy array 
    tree = cKDTree(centroids)
    return tree, grid_keys, centroids

#something is not working properly
#neighborhood of all grids 
#only compute the neighborhood for the particular grid
def neighborhood(cutoff:int, grid_number:int, tree, grid_keys:list, centroids):
    """
    This function creating a neighborhood.

    Args:
        cuttoff (int): The neighborhood cutoff distance in nm. 
        grid_number (int): The grid number for which to find neighborhoods.
        tree: The cKDTree object for spatial queries.
        grid_keys: List of grid numbers
        centroids: Numpy array of centroid coordinates for each grid 

    Returns:
        list: A list of nearby grids determined through cDKtree query.
    """
    #create neighborhood grid for each grid
    nearby_grids = tree.query_ball_point(centroids[grid_number], cutoff) #takes the grid point and pulls out the centroid coordinates
    #cutoff determines the points within the radius(cutoff), note this gives the indices into the array used to build the KDtree
    nearby_keys = [grid_keys[j] for j in nearby_grids] #turns the nearby_grids into a dictionary
    return nearby_keys






