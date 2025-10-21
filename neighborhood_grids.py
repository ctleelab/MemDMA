import numpy as np
from scipy.spatial import cKDTree



#something is not working properly
#neighborhood of all grids 
def neighborhood(grid_centroids, cutoff):
    grid_keys = list(grid_centroids.keys())  #creates a list of the keys of all centroids of the grid points
    centroids = np.array([grid_centroids[k][:3] for k in grid_keys]) #take the coordinates for each grid 
    tree = cKDTree(centroids)
    grid_neighborhood = {}
    #create neighborhood grid, but only if its not created already 
    for i, key in enumerate(grid_keys): #loop through all grid points
        nearby_grids = tree.query_ball_point(centroids[i], cutoff) #takes the grid point and pulls out the centroid coordinates
        #cutoff determines the points within the radius(cutoff), note this gives the indices into the array used to build the KDtree
        nearby_keys = [grid_keys[j] for j in nearby_grids if j != i]  # remove self
        grid_neighborhood[key] = nearby_keys
    return grid_neighborhood






