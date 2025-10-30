import numpy as np
from scipy.spatial import cKDTree


def cKDtree_grid_neighborhood(grid_centroids):
    grid_keys = list(grid_centroids.keys())  #creates a list of the keys of all centroids of the grid points
    centroids = np.array([grid_centroids[k][:3] for k in grid_keys]) #take the coordinates for each grid in grid_keys and creates a numpy array 
    tree = cKDTree(centroids)
    return tree, grid_keys, centroids

#something is not working properly
#neighborhood of all grids 
#only compute the neighborhood for the particular grid
def neighborhood(cutoff, grid_number, tree, grid_keys, centroids):

    grid_neighborhood = {}
    #create neighborhood grid for each grid
    nearby_grids = tree.query_ball_point(centroids[grid_number], cutoff) #takes the grid point and pulls out the centroid coordinates
    #cutoff determines the points within the radius(cutoff), note this gives the indices into the array used to build the KDtree
    nearby_keys = [grid_keys[j] for j in nearby_grids] #turns the nearby_grids into a dictionary
    grid_neighborhood = nearby_keys
    return grid_neighborhood






