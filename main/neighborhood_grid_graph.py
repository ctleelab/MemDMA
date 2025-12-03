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




def atom_mapping_centroids(triangle_centroids, triangle_atom_assignment):
    triangle_grid_all = []
    X_triangle_all = []
    Y_triangle_all = []
    Z_triangle_all = []
    count_for_triangle = []

    for key, value in triangle_centroids.items():
        centroid_num = key
        gl2_bead_coords = np.array(value)
        X = gl2_bead_coords[0]
        Y = gl2_bead_coords[1]
        Z = gl2_bead_coords[2]
        triangle_grid_all.append(centroid_num)
        X_triangle_all.append(X)
        Y_triangle_all.append(Y)
        Z_triangle_all.append(Z)
        count_for_triangle.append(triangle_atom_assignment.get(centroid_num, 0))  # Get count, default to 0 if not found

    # Create Plotly figure
    fig = go.Figure()

    # Add scatter plot for triangle centroids with count as hover text
    fig.add_trace(go.Scatter3d(
        x=X_triangle_all,
        y=Y_triangle_all,
        z=Z_triangle_all,
        mode='markers',
        name='Triangle Centroids',
        marker=dict(
            size=6,
            color=count_for_triangle,  # Color based on count
            colorscale='Viridis',  # Color scale
            opacity=0.8,
            colorbar=dict(title="Count")  # Add color bar for counts
        ),
        text=[f"Grid: {grid}, Count: {count}" for grid, count in zip(triangle_grid_all, count_for_triangle)],  # Hover text
        hoverinfo="text"
    ))

    fig.update_layout(
        title="Triangle Centroids with Atom Count",
        scene=dict(
            xaxis_title="X",
            yaxis_title="Y",
            zaxis_title="Z height",
        )
    )

    # Save as standalone HTML
    fig.write_html("triangle_centroids_with_counts.html")
    del fig
    print("Interactive plot with counts completed")


# goal graph gl2 bead and highlight the surrounding neighborhood
# input: gl2 beads
# input: grid neighborhoods 
def gl2_neighborhood(gl2_beads, triangle_centroids, gl2_grid_assignment):
    """
        graphs the neighborhood for each gl2 bead

    Args:
        gl2_beads (dict): A dictionary with GL2 (atom #): array([x y z]) coordinates of the GL2 bead
        gl2_grid_assignment (dict): A dictionary with GL2 (key): (closest grid_index, grid_coords of closest grid_index, neighborhood_list of the grid_index).
        triangle_centroids (dict): A dictionary with triangle number (key): [x y z] coordinates of the triangle centroid
    Returns:
        
    """
    # Extract triangle centroids
    triangle_grid_all = []
    X_triangle_all = []
    Y_triangle_all = []
    Z_triangle_all = []

    for key, value in triangle_centroids.items():
        centroid_num = key
        coords = np.array(value)
        X_triangle_all.append(coords[0])
        Y_triangle_all.append(coords[1])
        Z_triangle_all.append(coords[2])
        triangle_grid_all.append(centroid_num)

    # Extract GL2 bead data
    gl2_bead_coords = list(gl2_beads.values())[0]  # Take the first GL2 bead
    gl2_bead_x, gl2_bead_y, gl2_bead_z = gl2_bead_coords 

    # Get the grid assignment and neighborhood for the GL2 bead
    gl2_bead_key = list(gl2_beads.keys())[0]
    grid_assignment = gl2_grid_assignment[gl2_bead_key]
    neighborhood_list = grid_assignment[2]

    # Extract neighborhood centroids
    neighborhood_x = []
    neighborhood_y = []
    neighborhood_z = []

    for neighbor in neighborhood_list:
        neighbor_coords = triangle_centroids[neighbor]
        neighborhood_x.append(neighbor_coords[0])
        neighborhood_y.append(neighbor_coords[1])
        neighborhood_z.append(neighbor_coords[2])

    # Create Plotly figure
    fig = go.Figure()

    # Plot triangle centroids
    fig.add_trace(go.Scatter3d(
        x=X_triangle_all,
        y=Y_triangle_all,
        z=Z_triangle_all,
        mode='markers',
        name='Triangle Centroids',
        marker=dict(size=4, color='gray', opacity=0.6)
    ))

    # Plot GL2 bead
    fig.add_trace(go.Scatter3d(
        x=[gl2_bead_x],
        y=[gl2_bead_y],
        z=[gl2_bead_z],
        mode='markers',
        name='GL2 Bead',
        marker=dict(size=6, color='red', opacity=1.0)
    ))

    # Plot neighborhood centroids
    fig.add_trace(go.Scatter3d(
        x=neighborhood_x,
        y=neighborhood_y,
        z=neighborhood_z,
        mode='markers',
        name='Neighborhood',
        marker=dict(size=5, color='blue', opacity=0.8)
    ))

    fig.update_layout(
        title="GL2 Bead and Neighborhood",
        scene=dict(
            xaxis_title="X",
            yaxis_title="Y",
            zaxis_title="Z height",
        )
    )

    # Save as standalone HTML
    fig.write_html("gl2_neighborhood.html")
    del fig
    print("GL2 neighborhood plot completed")


#this graph should display the original atom coords and the projected atom coords
#input: atom coords, projected atom coords, a defined grid
#let's just try with gl2 beads
#def projected_atom_coords(gl2_beads, triangle_centroids, gl2_grid_assignment):

def projected_atom_grid(atom_coords_projections_absolute, neighborhood_centroid_coords_array, neighborhood_coords_array,atom_coords, neighborhood_norms_array,atom_num):
    """
        Graphs the projected atom coords for each grid. The coordinate for the centroid grid should match the projected atom coordination, by color.
        I also want to plot the triangle.

    Args:
        atom_coords_projections_absolute: List of projected atom coordinates.
        neighborhood_centroid_coords_array: List of neighborhood centroid coordinates.
        neighborhood_coords_array: List of vertices of (triangles) for each neighborhood.
        atom_cords: coordinates of atom 
    Returns:
        None
    """
    print("atom_coords:",atom_coords)
    # Prepare data for plotting
    X_projected = []
    Y_projected = []
    Z_projected = []
    X_centroid = []
    Y_centroid = []
    Z_centroid = []
    X_norms = []
    Y_norms =  []
    Z_norms = []
    vertex_X = []
    vertex_Y = []
    vertex_Z = []
    colors = []

    for i, atom in enumerate(atom_coords_projections_absolute):
        if i < 10 :
            # Projected atom coordinates
            X_projected.append(atom[0])
            Y_projected.append(atom[1])
            Z_projected.append(atom[2])

            # Centroid coordinates
            centroid = neighborhood_centroid_coords_array[i]
            X_centroid.append(centroid[0])
            Y_centroid.append(centroid[1])
            Z_centroid.append(centroid[2])

            # Triangle edges (3 vertices per triangle)
            grid = neighborhood_coords_array[i]
            for vertex in grid:
                vertex_X.append(vertex[0])
                vertex_Y.append(vertex[1])
                vertex_Z.append(vertex[2])

            #normal vectors of grids
            grid_norm = neighborhood_norms_array[i]
            X_norms.append(grid_norm[0])
            Y_norms.append(grid_norm[1])
            Z_norms.append(grid_norm[2])


            # Assign a unique color for each index
            colors.append(i)

    # Create Plotly figure
    fig = go.Figure()

    # Plot projected atom coordinates
    fig.add_trace(go.Scatter3d(
        x=X_projected,
        y=Y_projected,
        z=Z_projected,
        mode='markers',
        name='Projected Atoms',
        marker=dict(
            size=6,
            color=colors,  # Color based on index
            colorscale='Viridis',
            opacity=0.8,
            colorbar=dict(title="Index")
        )
    ))

    # Plot centroid coordinates
    fig.add_trace(go.Scatter3d(
        x=X_centroid,
        y=Y_centroid,
        z=Z_centroid,
        mode='markers',
        name='Centroid Grids',
        marker=dict(
            size=6,
            color=colors,  # Match color with projected atoms
            colorscale='Viridis',
            opacity=0.8
        )
    ))

    # Plot triangle edges
    fig.add_trace(go.Scatter3d(
        x=vertex_X,
        y=vertex_Y,
        z=vertex_Z,
        mode='lines',
        name='Triangle vertices',
        line=dict(color='blue', width=2)
    ))

    # # Plot normal vectors of grids
    # for i in range(len(X_centroid)):
    #     fig.add_trace(go.Scatter3d(
    #         x=[X_centroid[i], X_centroid[i] + X_norms[i]],
    #         y=[Y_centroid[i], Y_centroid[i] + Y_norms[i]],
    #         z=[Z_centroid[i], Z_centroid[i] + Z_norms[i]],
    #         mode='lines',
    #         name=f'Normal Vector {i}',
    #         line=dict(color='green', width=2)
    #     ))

    # Plot original atom 
    fig.add_trace(go.Scatter3d(
        x= [atom_coords[0]],
        y= [atom_coords[1]],
        z= [atom_coords[2]],
        mode='markers',
        name = 'Original atom coordinates',
        marker = dict(
            size = 6,
            color = 'red',
            opacity = .8
        )

    ))

    fig.update_layout(
        title="Projected Atom Coordinates and Centroid Grids",
        scene=dict(
            xaxis_title="X",
            yaxis_title="Y",
            zaxis_title="Z height",
        )
    )

    # Save as standalone HTML
    fig.write_html(f"projected_atom_grid_atom_num{atom_num}.html")
    print("Projected atom grid plot completed")






    


    