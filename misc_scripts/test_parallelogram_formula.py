import numpy as np
import matplotlib.pyplot as plt
import ctleelab_plothelper.plothelpers as ph


from matplotlib.patches import Circle

radius = 3
atom = np.array([2.6, 2.6])
vertices = np.array([[0, 0], [1, 0], [0, 1]])


fig, ax = ph.fixed_size_subplots(1, 1, subheight=3, subwidth=3)
ax.scatter(atom[0], atom[1], color="red")
circle = Circle(atom, radius=radius, color="red", alpha=0.3)
ax.add_patch(circle)


ax.plot(vertices[:, 0], vertices[:, 1], color="blue")
ax.set_aspect("equal")

vertex_to_atom_proj = atom - vertices

print(
    np.linalg.norm(vertex_to_atom_proj, axis=1),
    np.min(np.linalg.norm(vertex_to_atom_proj, axis=1)) <= radius,
)


edges = np.vstack(
    (
        vertices[1] - vertices[0],
        vertices[2] - vertices[1],
        vertices[0] - vertices[2],
    )
)

print(edges)
edge_lengths = np.linalg.norm(edges, axis=1)
print(edge_lengths)
edge_dots = np.einsum("ij,ij->i", vertex_to_atom_proj, edges)

print(edge_dots >= 0, edge_dots <= edge_lengths)

print(np.bitwise_and((edge_dots >= 0), (edge_dots <= edge_lengths)))

edge_projection_valid = np.bitwise_and((edge_dots >= 0), (edge_dots <= edge_lengths))
print(edge_projection_valid)

edge_valid = np.any(edge_projection_valid)


def cross2d(x, y):
    return x[..., 0] * y[..., 1] - x[..., 1] * y[..., 0]


dist_to_edge = cross2d(vertex_to_atom_proj, edges) / edge_lengths
print(dist_to_edge)

print((dist_to_edge <= radius) & edge_projection_valid)


print((dist_to_edge <= radius) & edge_valid)

# if np.any((dist_to_edge <= radius) & edge_projection_valid):
#     # print("edge valid")
#     increment_defect_count(leaflet, grid_idx, aliphatic)

plt.show()
