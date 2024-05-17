import os
import sys

import numpy as np

from utils import plot


# Ellipsoid: x^2 + y^2 = 1-2z^2
#   We only generate half of it for printing


def generate(
        width_mm=2,
        resolution_z_mm=0.1, resolution_xy_mm=0.1,
        range_x_cm=(-10, 10), range_y_cm=(0, 20), range_z_cm=(-10, 10)
        ):
    n_xy = int((range_x_cm[1]-range_x_cm[0])*np.sqrt(2) / resolution_xy_mm) + 1
    t = np.linspace(range_x_cm[0], range_x_cm[1], n_xy)

    x_coords = []
    y_coords = []
    z_coords = []
    # Parabolic Cylinder: x^2 = 2y
    #   For the extreme z values, parameterize (x, y) values
    for z in [range_z_cm[0], range_z_cm[1]]:
        x_coords.append(t)
        y_coords.append(t**2/2)
        z_coords.append(z*np.ones_like(t))

    X = np.concatenate(x_coords)
    Y = np.concatenate(y_coords)
    Z = np.concatenate(z_coords)

    nan_idcs = np.isnan(Y) | np.less(Y, range_y_cm[0]) | np.greater(Y, range_y_cm[1])
    X, Y, Z = X[~nan_idcs], Y[~nan_idcs], Z[~nan_idcs]

    # Normal direction at each point:
    normal = np.stack([X, -np.ones_like(Y), np.zeros_like(Z)], axis=-1)  # Normal: A x_0 + B -> [x_0, -1, 0]
    normal /= np.linalg.norm(normal, axis=-1)[..., np.newaxis]

    # Generate polyhedrons for OpenSCAD
    #   Each point is a list of 3 coordinates
    points_outer = np.stack([X, Y, Z], axis=-1)
    points_inner = points_outer - normal * (width_mm / 10)
    points_cm = np.concatenate([points_outer, points_inner], axis=0)
    points_mm = points_cm * 10
    plot(points_outer, points_inner)

    #   Each face is a collection of 6 points:
    #       bot_x_first, bot_x_second, bot_y_second, bot_y_first
    #       top_x_first, top_x_second, top_y_second, top_y_first
    n_xy = X.size // 2
    idcs = np.arange(n_xy)
    faces_side_outer = np.stack([idcs[1:], idcs[:-1], idcs[:-1] + n_xy, idcs[1:] + n_xy], axis=-1)
    faces_sides = np.concatenate([faces_side_outer, faces_side_outer + n_xy*2], axis=0)

    faces_caps = np.concatenate([
        np.stack([idcs[1:], idcs[:-1], idcs[:-1] + n_xy*2, idcs[1:] + n_xy*2], axis=-1),
        np.stack([idcs[1:]+n_xy, idcs[:-1]+n_xy, idcs[:-1]+n_xy + n_xy*2, idcs[1:]+n_xy + n_xy*2], axis=-1),
        np.asarray([[0, n_xy, 3 * n_xy, 2 * n_xy], [n_xy - 1, 2 * n_xy - 1, 4 * n_xy - 1, 3 * n_xy - 1]])
    ])

    # Check points are within the range
    assert np.all((range_x_cm[0] <= points_cm[:, 0]) & (points_cm[:, 0] <= range_x_cm[1])), 'Out of bounds'
    assert np.all((range_y_cm[0] <= points_cm[:, 1]) & (points_cm[:, 1] <= range_y_cm[1])), 'Out of bounds'
    assert np.all((range_z_cm[0] <= points_cm[:, 2]) & (points_cm[:, 2] <= range_z_cm[1])), 'Out of bounds'

    np.set_printoptions(threshold=sys.maxsize)
    os.makedirs(f"scad_{width_mm}mm", exist_ok=True)
    with open(f"scad_{width_mm}mm/paraboliccylinder.scad", "w") as f:
        f.write(f"points = {np.array2string(points_mm, separator=', ')};\n\n")
        f.write(f"faces_sides = {np.array2string(faces_sides, separator=', ')};\n\n")
        f.write(f"faces_caps = {np.array2string(faces_caps, separator=', ')};\n\n")
        f.write("polyhedron(points, concat(faces_sides, faces_caps), convexity=1);")


def radius(z):
    """ Radius of the conic section at height z. """
    return np.sqrt(9 - z**2)


if __name__ == "__main__":
    generate()