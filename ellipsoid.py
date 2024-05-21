import os
import sys

import numpy as np

from utils import plot


def generate(
        a=1, b=2, c=1, r=5,
        width_mm=2,
        resolution_z_mm=0.1, resolution_xy_mm=0.1,
        range_x_cm=(-10, 10), range_y_cm=(-10, 10), range_z_cm=(0, 10),
        ):
    """ Ellipsoid: a*x^2 + b*y^2 +c*z^2 = r^2.
    We only generate half of it for printing
    """
    min_z = max([range_z_cm[0], range_x_cm[0] * np.sqrt(a), range_y_cm[0] * np.sqrt(b), -r/np.sqrt(c)])
    max_z = min([range_z_cm[1], range_x_cm[1] * np.sqrt(a), range_y_cm[1] * np.sqrt(b), +r/np.sqrt(c)])
    max_n_z = int((max_z - min_z) / resolution_z_mm) + 1
    max_perimeter_cm = 2 * np.pi * np.nanmax(np.sqrt(1 - 2*np.linspace(range_z_cm[0], range_z_cm[1], max_n_z)**2))
    n_xy = int(max_perimeter_cm / resolution_xy_mm) + 1
    theta = np.linspace(0, 2 * np.pi, n_xy)

    x_coords = []
    y_coords = []
    z_coords = []
    #   For fixed z values, parameterize (x, y) values
    n_z = 0
    for z in np.linspace(min_z, max_z, max_n_z):
        r_z = np.sqrt(r*r - c*z*z)
        if not np.isnan(r):
            n_z += 1
            x_coords.append(r_z*np.cos(theta)/np.sqrt(a))
            y_coords.append(r_z*np.sin(theta)/np.sqrt(b))
            z_coords.append(z*np.ones_like(theta))

    x_coords.append([0])
    y_coords.append([0])
    z_coords.append([+r/np.sqrt(c)])

    X = np.concatenate(x_coords)
    Y = np.concatenate(y_coords)
    Z = np.concatenate(z_coords)

    # Normal direction at each point:
    normal = np.stack([a*X, b*Y, 2*Z], axis=-1)  # Normal: A x_0 + B -> [x_0, y_0, 2*z_0]
    normal /= np.linalg.norm(normal, axis=-1)[..., np.newaxis]

    # Generate polyhedrons for OpenSCAD
    #   Each point is a list of 3 coordinates
    points_outer = np.stack([X, Y, Z], axis=-1)
    points_inner = points_outer - normal * (width_mm / 10)
    # Correct the range of z since normal vector can have contribution in this direction
    points_inner[:, 2] = np.maximum(points_inner[:, 2], range_z_cm[0])
    points_inner[:, 2] = np.minimum(points_inner[:, 2], range_z_cm[1])
    points_cm = np.concatenate([points_outer, points_inner], axis=0)
    points_mm = 10*points_cm
    plot(points_outer, points_inner)

    #   Each face is a collection of 6 points:
    #       bot_x_first, bot_x_second, bot_y_second, bot_y_first
    #       top_x_first, top_x_second, top_y_second, top_y_first
    idcs = np.arange(n_xy)
    faces_first_layer = np.concatenate([
        np.stack([idcs[1:], idcs[:-1], idcs[:-1] + n_xy, idcs[1:] + n_xy], axis=-1),
        np.stack([idcs[-1:], idcs[:1], idcs[:1] + n_xy, idcs[-1:] + n_xy], axis=-1)
    ])
    faces_side_outer = np.concatenate([faces_first_layer + n_xy * i for i in range(n_z-1)], axis=0)
    faces_sides = np.concatenate([faces_side_outer, faces_side_outer + n_xy*n_z+1], axis=0)

    idx_last = n_xy * n_z
    faces_cap_outer = np.concatenate([
        np.stack([idcs[1:]+n_xy*(n_z-1), idcs[:-1]+n_xy*(n_z-1), idx_last*np.ones_like(idcs[1:])], axis=-1),
        np.stack([idcs[-1:]+n_xy*(n_z-1), idcs[:1]+n_xy*(n_z-1), idx_last*np.ones_like(idcs[-1:])], axis=-1)
    ])
    faces_caps = np.concatenate([faces_cap_outer, faces_cap_outer + n_xy*n_z+1], axis=0)

    faces_bottom = np.stack([idcs[1:], idcs[:-1], idcs[:-1] + n_xy*n_z+1, idcs[1:] + n_xy*n_z+1], axis=-1)

    # Check points are within the range
    assert np.all((range_x_cm[0] <= points_cm[:, 0]) & (points_cm[:, 0] <= range_x_cm[1])), 'Out of bounds'
    assert np.all((range_y_cm[0] <= points_cm[:, 1]) & (points_cm[:, 1] <= range_y_cm[1])), 'Out of bounds'
    assert np.all((range_z_cm[0] <= points_cm[:, 2]) & (points_cm[:, 2] <= range_z_cm[1])), 'Out of bounds'

    np.set_printoptions(threshold=sys.maxsize)
    os.makedirs(f"scad_{width_mm:.01f}mm", exist_ok=True)
    with open(f"scad_{width_mm:.01f}mm/ellipsoid.scad", "w") as f:
        f.write(f"points = {np.array2string(points_mm, separator=', ')};\n\n")
        f.write(f"faces_sides = {np.array2string(faces_sides, separator=', ')};\n\n")
        f.write(f"faces_caps = {np.array2string(faces_caps, separator=', ')};\n\n")
        f.write(f"faces_bottom = {np.array2string(faces_bottom, separator=', ')};\n\n")
        f.write("polyhedron(points, concat(faces_sides, faces_caps, faces_bottom), convexity=1);")



if __name__ == "__main__":
    generate()