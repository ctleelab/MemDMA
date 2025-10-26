#general plot of defects

import matplotlib.pyplot as plt

import plotly.graph_objects as go

import pickle

def bar_graph_defects(pickle_file):
    with open(pickle_file, "rb") as f:
        neighborhood_grid_id_upper = pickle.load(f)

    below_1 = 0
    above_1 = 0
    total_grids = len(neighborhood_grid_id_upper)
    for key, value in neighborhood_grid_id_upper.items():
        sum_atoms = value
        if sum_atoms >= 1:
            above_1 += 1
        if sum_atoms < 1:
            below_1 += 1
    #normalize by total number of grids
    below_1_norm = below_1/total_grids
    above_1_norm = above_1/total_grids
    
    # Create Plotly bar chart
    fig = go.Figure(
        data=[
            go.Bar(
                x=["defect", "no defect"],
                y=[below_1_norm, above_1_norm],
                marker_color=["#636EFA", "#EF553B"],
            )
        ]
    )

    fig.update_layout(
        title="Distribution of Grid Assignments",
        xaxis_title="Atom Mapping Sum Category",
        yaxis_title="Normalized count of Triangles",
        template="plotly_white",
    )

    # Show plot in browser or save to HTML
    fig.write_html("DPPC_triangle_atom_assignment_upper_strain2.html")
    print("✅ Bar graph saved'")

    return below_1_norm, above_1_norm
    


below_1_norm_DOPC_flat, above_1_norm_DOPC_flat = bar_graph_defects("DOPC_triangle_atom_assignment_upper_flat.pickle")

below_1_norm_DOPC_strain_2, above_1_norm_DOPC_strain_2 = bar_graph_defects("DOPC_triangle_atom_assignment_upper_buckled.pickle")

below_1_norm_DPPC_flat, above_1_norm_DPPC_flat = bar_graph_defects("DPPC_triangle_atom_assignment_upper_flat.pickle")

below_1_norm_DPPC_strain_2, above_1_norm_DPPC_strain_2 = bar_graph_defects("DPPC_triangle_atom_assignment_upper_strain2.pickle")




def bar_graph_defects_all_simulations(below_1_norm_DOPC_flat, below_1_norm_DOPC_strain_2, below_1_norm_DPPC_flat, below_1_norm_DPPC_strain_2):
    
    # Create Plotly bar chart
    fig = go.Figure(
        data=[
            go.Bar(
                x=["DOPC flat","DPPC flat","DOPC strain .2","DPPC strain .2"],
                y=[below_1_norm_DOPC_flat,below_1_norm_DPPC_flat,below_1_norm_DOPC_strain_2, below_1_norm_DPPC_strain_2],
            )
        ]
    )

    fig.update_layout(
        title="Defects comparison across simulations",
        xaxis_title="simulation type",
        yaxis_title="Percentage of grids that have defects",
        template="plotly_white",
    )

    # Show plot in browser or save to HTML
    fig.write_html("summary_simulations_defects.html")
    print("✅ Bar graph saved'")

bar_graph_defects_all_simulations(below_1_norm_DOPC_flat,below_1_norm_DOPC_strain_2,below_1_norm_DPPC_flat,below_1_norm_DPPC_strain_2)
