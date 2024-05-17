import os
import sys

import numpy as np

from utils import plot


def generate(
        a=1, b=2,
        width_mm=2,
        resolution_z_mm=0.2, resolution_xy_mm=0.2,
        range_x_cm=(-10, 10), range_y_cm=(-10, 10), range_z_cm=(-10, 10)
        ):
    """ Hyperboloid a*x^2 + b*y^2 - z^2 = 2*(3*3). """
    min_z = max([range_z_cm[0], -np.sqrt(range_x_cm[0]**2 * a - 18), -np.sqrt(range_y_cm[0]**2 * b - 18)])
    max_z = min([range_z_cm[1], np.sqrt(range_x_cm[1]**2 * a - 18), np.sqrt(range_y_cm[1]**2 * b - 18)])
    n_z = int((max_z - min_z) / resolution_z_mm) + 1
    max_perimeter_cm = 2 * np.pi * np.max(radius(np.linspace(min_z, max_z, n_z)))
    n_xy = int(max_perimeter_cm / resolution_xy_mm) + 1
    theta = np.linspace(0, 2 * np.pi, n_xy)

    x_coords = []
    y_coords = []
    z_coords = []
    #   For fixed z values, parameterize (x, y) values
    for z in np.linspace(min_z, max_z, n_z):
        r = radius(z)
        x_coords.append(r*np.cos(theta)/np.sqrt(a))
        y_coords.append(r*np.sin(theta)/np.sqrt(b))
        z_coords.append(z*np.ones_like(theta))

    X = np.concatenate(x_coords)
    Y = np.concatenate(y_coords)
    Z = np.concatenate(z_coords)

    # Normal direction at each point:
    normal = np.stack([a*X, b*Y, -np.ones_like(Z)], axis=-1)   # Normal: A x_0 + B -> [a*x_0, b*y_0, -1]
    normal /= np.linalg.norm(normal, axis=-1)[..., np.newaxis]

    outer = np.stack([X, Y, Z], axis=-1)
    inner = outer - normal * (width_mm / 10)

    plot(outer, inner)

    # Generate polyhedrons for OpenSCAD
    #   Each point is a list of 3 coordinates
    points_cm = np.concatenate([outer, inner], axis=0)
    # Correct the ranges since normal vector can have contributions in each direction
    points_cm[:, 0] = np.maximum(points_cm[:, 0], range_x_cm[0])
    points_cm[:, 0] = np.minimum(points_cm[:, 0], range_x_cm[1])
    points_cm[:, 1] = np.maximum(points_cm[:, 1], range_y_cm[0])
    points_cm[:, 1] = np.minimum(points_cm[:, 1], range_y_cm[1])
    points_cm[:, 2] = np.maximum(points_cm[:, 2], range_z_cm[0])
    points_cm[:, 2] = np.minimum(points_cm[:, 2], range_z_cm[1])
    points_mm = points_cm * 10

    #   Each face is a collection of 6 points:
    #       bot_x_first, bot_x_second, bot_y_second, bot_y_first
    #       top_x_first, top_x_second, top_y_second, top_y_first
    idcs = np.arange(n_xy)
    faces_first_layer = np.concatenate([
        np.stack([idcs[1:], idcs[:-1], idcs[:-1] + n_xy, idcs[1:] + n_xy], axis=-1),
        np.stack([idcs[-1:], idcs[:1], idcs[:1] + n_xy, idcs[-1:] + n_xy], axis=-1)
    ])
    faces_outer = np.concatenate([faces_first_layer + n_xy * i for i in range(n_z-1)], axis=0)
    face_inner = faces_outer + n_xy * n_z

    faces_join_bottom = np.concatenate([
        np.stack([idcs[1:], idcs[:-1], idcs[:-1] + n_xy * n_z, idcs[1:] + n_xy * n_z], axis=-1),
        np.stack([idcs[-1:], idcs[:1], idcs[:1] + n_xy * n_z, idcs[-1:] + n_xy * n_z], axis=-1)
    ])
    faces_join_top = faces_join_bottom + n_xy * (n_z - 1)

    faces = np.concatenate([faces_outer, face_inner, faces_join_bottom, faces_join_top], axis=0)

    # Check points are within the range
    assert np.all((range_x_cm[0] <= points_cm[:, 0]) & (points_cm[:, 0] <= range_x_cm[1])), 'Out of bounds'
    assert np.all((range_y_cm[0] <= points_cm[:, 1]) & (points_cm[:, 1] <= range_y_cm[1])), 'Out of bounds'
    assert np.all((range_z_cm[0] <= points_cm[:, 2]) & (points_cm[:, 2] <= range_z_cm[1])), 'Out of bounds'

    np.set_printoptions(threshold=sys.maxsize)
    os.makedirs(f"scad_{width_mm:.01f}mm", exist_ok=True)
    with open(f"scad_{width_mm:.01f}mm/hyperboloid.scad", "w") as f:
        f.write(f"points = {np.array2string(points_mm, separator=', ')};\n\n")
        f.write(f"faces = {np.array2string(faces, separator=', ')};\n\n")
        f.write("polyhedron(points, faces, convexity=1);")


def radius(z):
    """ Radius of the conic section at height z. """
    return np.sqrt(18 + z**2)


if __name__ == "__main__":
    generate()