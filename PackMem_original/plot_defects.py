import numpy as np
import plotly.graph_objects as go
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from pathlib import Path

testing_path = Path("/home/casakurai/scratch/packmem/PackMem_original/testing")

filename = "DOPC_flat_DefectGrid_Upper.txt"



# def plot_v1(filename): 
#     grid_spacing = 1.0  # Å (PackMem default)
#     xs, ys, areas = [], [], []


# def plot_v1(filename): 
#     grid_spacing = 1.0  # Å (PackMem default)
#     xs, ys, areas = [], [], []

#     with open(filename) as f:
#         for line in f:
#             if not line.strip() or line.startswith("#"):
#                 continue
#             parts = line.split()
#             if len(parts) != 4:
#                 continue
#             try:
#                 area = float(parts[1])
#                 x_grid = float(parts[2])
#                 y_grid = float(parts[3])
#             except ValueError:
#                 continue
#     with open(filename) as f:
#         for line in f:
#             if not line.strip() or line.startswith("#"):
#                 continue
#             parts = line.split()
#             if len(parts) != 4:
#                 continue
#             try:
#                 area = float(parts[1])
#                 x_grid = float(parts[2])
#                 y_grid = float(parts[3])
#             except ValueError:
#                 continue

#             xs.append(x_grid * grid_spacing)
#             ys.append(y_grid * grid_spacing)
#             areas.append(area)
#             xs.append(x_grid * grid_spacing)
#             ys.append(y_grid * grid_spacing)
#             areas.append(area)

#     xs = np.array(xs)
#     ys = np.array(ys)
#     areas = np.array(areas)
#     xs = np.array(xs)
#     ys = np.array(ys)
#     areas = np.array(areas)

#     # scale marker sizes
#     marker_sizes = areas * 4
#     # scale marker sizes
#     marker_sizes = areas * 4

#     fig = go.Figure()
#     fig = go.Figure()

#     fig.add_trace(
#         go.Scatter(
#             x=xs,
#             y=ys,
#             mode="markers",
#             marker=dict(
#                 size=marker_sizes,
#                 color=areas,
#                 colorscale="Viridis",
#                 sizemode="area",
#                 showscale=True,
#                 colorbar=dict(title="Defect area (grid cells)"),
#             ),
#             text=[
#                 f"Area: {a:.0f} Å²<br>x={x:.1f} Å, y={y:.1f} Å"
#                 for a, x, y in zip(areas, xs, ys)
#             ],
#             hoverinfo="text",
#         )
#     )
#     fig.add_trace(
#         go.Scatter(
#             x=xs,
#             y=ys,
#             mode="markers",
#             marker=dict(
#                 size=marker_sizes,
#                 color=areas,
#                 colorscale="Viridis",
#                 sizemode="area",
#                 showscale=True,
#                 colorbar=dict(title="Defect area (grid cells)"),
#             ),
#             text=[
#                 f"Area: {a:.0f} Å²<br>x={x:.1f} Å, y={y:.1f} Å"
#                 for a, x, y in zip(areas, xs, ys)
#             ],
#             hoverinfo="text",
#         )
#     )

#     fig.update_layout(
#         title="Packing Defects in Membrane Plane (Å Coordinates)",
#         xaxis_title="X (Å)",
#         yaxis_title="Y (Å)",
#         template="plotly_white",
#     )
#     fig.update_layout(
#         title="Packing Defects in Membrane Plane (Å Coordinates)",
#         xaxis_title="X (Å)",
#         yaxis_title="Y (Å)",
#         template="plotly_white",
#     )

#     fig.write_html(f"{testing_path}/packmem_defects_xy_angstrom.html")
#     print("Saved packmem_defects_xy_angstrom.html")


# #plot per grid defect 


# def plot_binary_defects_red0(filetxt, leaflet_name="Up", shift_coords=True):
#     """
#     Plot binarized packing defects where 0 is red and 1 is white.
    
#     Parameters:
#         filetxt (str): path to the .txt file (x, y, leaflet, defect)
#         leaflet_name (str): 'Up' or 'Lo', used for title
#         shift_coords (bool): shift x,y to be positive
#     """
#     # Load data
#     df = pd.read_csv(filetxt, sep="\t")
    
#     # Binarize defect column: 1 if defect >0 else 0
#     df['defect'] = df['defect'].apply(lambda x: 1 if x > 0 else 0)
    
#     # Optional: shift coordinates to positive
#     if shift_coords:
#         df['x'] = df['x'] - df['x'].min()
#         df['y'] = df['y'] - df['y'].min()
    
#     # Pivot to 2D grid
#     grid = df.pivot(index='y', columns='x', values='defect')
    
#     # Create custom colormap: 0 -> red, 1 -> white
#     cmap = mcolors.ListedColormap(['red', 'white'])
    
#     plt.figure(figsize=(8,6))
#     plt.imshow(grid, origin='lower', cmap=cmap, interpolation='none',
#                extent=[df['x'].min(), df['x'].max(), df['y'].min(), df['y'].max()])
    
#     # Add grid legend manually
#     import matplotlib.patches as mpatches
#     red_patch = mpatches.Patch(color='red', label='No defect (0)')
#     white_patch = mpatches.Patch(color='white', label='Defect (1)')
#     plt.legend(handles=[red_patch, white_patch], loc='upper right')
    
#     plt.xlabel('X')
#     plt.ylabel('Y')
#     plt.title(f'{leaflet_name} Leaflet Packing Defects (0=red, 1=white)')
#     plt.savefig(f"{filetxt}_defects.png")
#     plt.show()


#scatterplot of defects, where color represents defect type
def plot1(txtfile):
    # read the text file
    df = pd.read_csv(
        txtfile,
        delim_whitespace=True
    )

    # map colors: red if Value <= 1 white otherwise
    colors = np.where(df["Value"] <= 1, "red", "white")

    plt.figure(figsize=(5, 5))
    plt.scatter(
    df["X"],
    df["Y"],
    c=colors)

    plt.xlabel("X")
    plt.ylabel("Y")
    plt.savefig("Packing-defects.png")
    plt.close()

plot1(f"{testing_path}/{filename}")



