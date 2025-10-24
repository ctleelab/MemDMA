import plotly.graph_objects as go
import numpy as np


#input upper_neighborhood_cKDtree.csv
#key is the grid
#list is a list of grid

#what grids are grouping together in space..? 
#input also the centroid 




#can I plot the GL2 bead in space? based on the grid it assigns to..?
#then can I look at which grids are being assigned to it

#take GL2 bead, grab the GL2 closest grid, grab the coords, then plot coords in 3D space

#{gl2 atom_number: grid_number, (coords of grid), list of neighborhood grids}

def plot_gl2_grid_assignment(gl2_grid_assignment_upper, gl2_grid_assignment_lower, gl2_bead): 


    grid_num_upper_all = []
    X_upper_all = [] 
    Y_upper_all = []
    Z_upper_all = []

    for key, value in gl2_grid_assignment_upper.items():
        grid_upper_number = value[0]
        grid_upper_coords = value[1]
        #print("gl2 bead:", key, "grid number:", grid_upper_number, "grid_coords:", grid_upper_coords) #all gl2 beads are mapping to the same grid..
        X_up = grid_upper_coords[0]
        Y_up = grid_upper_coords[1]
        Z_up = grid_upper_coords[2]

        grid_num_upper_all.append(grid_upper_number)
        X_upper_all.append(X_up)
        Y_upper_all.append(Y_up)
        Z_upper_all.append(Z_up)



    grid_num_lower_all = []
    X_lower_all = [] 
    Y_lower_all = []
    Z_lower_all = []

    for key2, value2 in gl2_grid_assignment_lower.items():
        grid_lower_number = value2[0]
        grid_lower_coords = value2[1]
        #something is wrong with the grid_coords. They are all almost alike.
        X_low = grid_lower_coords[0]
        Y_low = grid_lower_coords[1]
        Z_low = grid_lower_coords[2]
        grid_num_lower_all.append(grid_lower_number)
        X_lower_all.append(X_low)
        Y_lower_all.append(Y_low)
        Z_lower_all.append(Z_low)


    gl2_bead_all = []
    X_gl2_all = []
    Y_gl2_all = []
    Z_gl2_all = []

    for key1, value1 in gl2_bead.items():
        gl2_bead_num = key1
        
        gl2_bead_coords = np.array(value1)
        #if I multiply by 10, it brings it to A units 
        X_gl2 = gl2_bead_coords[0]
        Y_gl2 = gl2_bead_coords[1]
        Z_gl2 = gl2_bead_coords[2]
        gl2_bead_all.append(gl2_bead_num)
        X_gl2_all.append(X_gl2)
        Y_gl2_all.append(Y_gl2)
        Z_gl2_all.append(Z_gl2)



    # Create Plotly figure
    fig = go.Figure()

    #set min and max height (z-cord) between the leaflets and use that to set the color scale 


    # Create Plotly figure
    # Black is the coordinates of the grids
    # Red is the coordiantes of the beads 
    fig.add_trace(go.Scatter3d(
        x = X_upper_all,
        y = Y_upper_all,
        z = Z_upper_all, 
        mode = 'markers',
        name = 'coords of grids upper',
        marker=dict(size=4, color = "black", opacity=0.8)
    ))

    fig.add_trace(go.Scatter3d(
        x = X_lower_all,
        y = Y_lower_all,
        z = Z_lower_all, 
        mode = 'markers',
        name = 'coords of grids lower',
        marker=dict(size=4, color = "blue", opacity=0.8)
    ))


    fig.add_trace(go.Scatter3d(
        x = X_gl2_all,
        y = Y_gl2_all,
        z = Z_gl2_all, 
        mode = 'markers',
        name = 'coords of GL2 beads',
        marker=dict(size=4, color = "red", opacity=0.8)
    ))

    fig.update_layout(
        title=f"Membrane Leaflet upper",
        scene=dict(
            xaxis_title="X",
            yaxis_title="Y",
            zaxis_title="Z height",
        )
    )

    
    # Save as standalone HTML
    fig.write_html("testing.html")
    del fig
    print("plot completed")




def plot_triangle_grid_assignment(triangle_centroids_upper, triangle_centroids_lower): 

    triangle_grid_upper_all = []
    X_triangle_upper_all = []
    Y_triangle_upper_all = []
    Z_triangle_upper_all = []

    for key, value in triangle_centroids_upper.items():
        centroid_num_upper = key
    
        gl2_bead_coords = np.array(value)
        X = gl2_bead_coords[0]
        Y = gl2_bead_coords[1]
        Z = gl2_bead_coords[2]
        triangle_grid_upper_all.append(centroid_num_upper)
        X_triangle_upper_all .append(X)
        Y_triangle_upper_all .append(Y)
        Z_triangle_upper_all.append(Z)


    triangle_grid_lower_all = []
    X_triangle_lower_all = []
    Y_triangle_lower_all = []
    Z_triangle_lower_all = []

    for key1, value1 in triangle_centroids_lower.items():
        centroid_num_lower = key1
    
        gl2_bead_coords = np.array(value1)
        X2 = gl2_bead_coords[0]
        Y2 = gl2_bead_coords[1]
        Z2 = gl2_bead_coords[2]
        triangle_grid_lower_all .append(centroid_num_lower)
        X_triangle_lower_all.append(X2)
        Y_triangle_lower_all.append(Y2)
        Z_triangle_lower_all.append(Z2)



    # Create Plotly figure
    fig = go.Figure()

    #set min and max height (z-cord) between the leaflets and use that to set the color scale 


    # Create Plotly figure
    # Black is the coordinates of the grids
    # Red is the coordiantes of the beads 
    fig.add_trace(go.Scatter3d(
        x = X_triangle_upper_all,
        y = Y_triangle_upper_all,
        z = Z_triangle_upper_all, 
        mode = 'markers',
        name = 'coords triangle centroids upper',
        marker=dict(size=4, color = "black", opacity=0.8)
    ))

    fig.add_trace(go.Scatter3d(
        x = X_triangle_lower_all,
        y = Y_triangle_lower_all,
        z = Z_triangle_lower_all, 
        mode = 'markers',
        name = 'coords triangle centroids lower',
        marker=dict(size=4, color = "red", opacity=0.8)
    ))

    fig.update_layout(
        title=f"Membrane triangle coords",
        scene=dict(
            xaxis_title="X",
            yaxis_title="Y",
            zaxis_title="Z height",
        )
    )

    
    # Save as standalone HTML
    fig.write_html("troubleshooting_triangle_centroids.html")
    del fig
    print("plot completed")







 



 
