import os
import sys

import numpy as np

from utils import plot


# Cone: z^2 = a*x^2 + b*y^2
#   We only generate half of it for printing

width_mm = 2
resolution_z_mm = 0.1
resolution_xy_mm = 0.1

range_x_cm = [-10, 10]
range_y_cm = [-10, 10]
range_z_cm = [-10, 0]


def generate(a=1, b=1):
    n_z = int((range_z_cm[1] - range_z_cm[0]) / resolution_z_mm) + 1
    max_perimeter_cm = 2 * np.pi * -range_z_cm[0]
    n_xy = int(max_perimeter_cm / resolution_xy_mm) + 1
    theta = np.linspace(0, 2 * np.pi, n_xy)

    # Cone: z^2 = a*x^2 + b*y^2
    #   For fixed z values, parameterize (x, y) values
    x_coords = [range_z_cm[0]*np.cos(theta)/np.sqrt(a)] + [[0]]
    y_coords = [range_z_cm[0]*np.sin(theta)/np.sqrt(b)] + [[0]]
    z_coords = [range_z_cm[0]*np.ones_like(theta)] + [[0]]

    X = np.concatenate(x_coords)
    Y = np.concatenate(y_coords)
    Z = np.concatenate(z_coords)

    # Normal direction at each point:
    normal = np.stack([a*X, b*Y, -Z], axis=-1)  # Normal: A x_0 + B -> [a*x_0, b*y_0, -z_0]
    normal /= np.linalg.norm(normal, axis=-1)[..., np.newaxis]
    normal[-1, :] = [0, 0, 1]   # Last point is the vertex at the tip of the cone

    # Generate polyhedrons for OpenSCAD
    #   Each point is a list of 3 coordinates
    points_outer = np.stack([X, Y, Z], axis=-1)
    points_inner = points_outer - normal * (width_mm / 10)
    # Correct the range of z since normal vector can have contribution in this direction
    points_inner[:, 2] = np.maximum(points_inner[:, 2], range_z_cm[0])
    points_inner[:, 2] = np.minimum(points_inner[:, 2], range_z_cm[1])
    points = np.concatenate([points_outer, points_inner], axis=0)
    plot(points_outer, points_inner)

    #   Each face is a collection of 6 points:
    #       bot_x_first, bot_x_second, bot_y_second, bot_y_first
    #       top_x_first, top_x_second, top_y_second, top_y_first
    idcs = np.arange(n_xy)

    idx_last = n_xy
    faces_cap_outer = np.concatenate([
        np.stack([idcs[1:], idcs[:-1], idx_last*np.ones_like(idcs[1:])], axis=-1),
        np.stack([idcs[-1:], idcs[:1], idx_last*np.ones_like(idcs[-1:])], axis=-1)
    ])
    faces_caps = np.concatenate([faces_cap_outer, faces_cap_outer + n_xy+1], axis=0)

    faces_bottom = np.stack([idcs[1:], idcs[:-1], idcs[:-1] + n_xy+1, idcs[1:] + n_xy+1], axis=-1)

    np.set_printoptions(threshold=sys.maxsize)
    os.makedirs("scad", exist_ok=True)
    with open("scad/cone_{a:.02f}_{b:.02f}.scad", "w") as f:
        f.write(f"points = {np.array2string(points, separator=', ')};\n\n")
        f.write(f"faces_caps = {np.array2string(faces_caps, separator=', ')};\n\n")
        f.write(f"faces_bottom = {np.array2string(faces_bottom, separator=', ')};\n\n")
        f.write("polyhedron(points, concat(faces_caps, faces_bottom), convexity=1);")



if __name__ == "__main__":
    generate(0.5, 0.5)
    generate(0.5, 1)
    generate(2, 1)