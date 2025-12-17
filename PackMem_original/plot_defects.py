import numpy as np
import plotly.graph_objects as go

filename = "Total__deep.txt"

grid_spacing = 1.0  # Å (PackMem default)

xs, ys, areas = [], [], []

with open(filename) as f:
    for line in f:
        if not line.strip() or line.startswith("#"):
            continue
        parts = line.split()
        if len(parts) != 4:
            continue
        try:
            area = float(parts[1])
            x_grid = float(parts[2])
            y_grid = float(parts[3])
        except ValueError:
            continue

        xs.append(x_grid * grid_spacing)
        ys.append(y_grid * grid_spacing)
        areas.append(area)

xs = np.array(xs)
ys = np.array(ys)
areas = np.array(areas)

# scale marker sizes
marker_sizes = areas * 4

fig = go.Figure()

fig.add_trace(
    go.Scatter(
        x=xs,
        y=ys,
        mode="markers",
        marker=dict(
            size=marker_sizes,
            color=areas,
            colorscale="Viridis",
            sizemode="area",
            showscale=True,
            colorbar=dict(title="Defect area (grid cells)"),
        ),
        text=[
            f"Area: {a:.0f} Å²<br>x={x:.1f} Å, y={y:.1f} Å"
            for a, x, y in zip(areas, xs, ys)
        ],
        hoverinfo="text",
    )
)

fig.update_layout(
    title="Packing Defects in Membrane Plane (Å Coordinates)",
    xaxis_title="X (Å)",
    yaxis_title="Y (Å)",
    template="plotly_white",
)

fig.write_html("packmem_defects_xy_angstrom.html")
print("Saved packmem_defects_xy_angstrom.html")
