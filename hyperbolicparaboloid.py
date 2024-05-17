import os
import sys

import numpy as np

from utils import plot


def generate(
        width_mm=2,
        resolution_z_mm=0.2, resolution_xy_mm=0.2,
        range_x_cm=(-10, 10), range_y_cm=(-10, 10), range_z_cm=(-10, 10)
        ):
    n_z = int((range_z_cm[1] - 0) / resolution_z_mm) + 1
    n_xy = int((range_x_cm[1]-range_x_cm[0]) * np.sqrt(2) / resolution_xy_mm)//2 * 2 + 1
    theta = np.linspace(0, 2 * np.pi, n_xy)

    x_coords = []
    y_coords = []
    z_coords = []
    # Hyperbolic paraboloid: x^2 - y^2 = 2*z
    #   For fixed z values, parameterize (x, y) values
    for z in np.linspace(0, range_z_cm[1], n_z):
        if z == 0:
            t = np.concatenate([
                -np.linspace((-min(range_x_cm[0], range_y_cm[0]))**0.5, 0, (n_xy-1)//2, endpoint=False)**2,
                np.asarray([0]),
                np.linspace((max(range_x_cm[1], range_y_cm[1]))**0.5, 0, (n_xy-1)//2, endpoint=False)[::-1]**2
            ])

            x_coords.append(np.abs(t))
            y_coords.append(t)
            z_coords.append(z*np.ones_like(t))
        else:
            theta_max = np.maximum(
                np.arccosh(range_x_cm[1] / np.sqrt(2*z)),
                np.arccosh(range_y_cm[1] / np.sqrt(2 * z))
            )
            theta_min = np.arcsinh(range_y_cm[0] / np.sqrt(2*z))
            theta = np.linspace(theta_min, theta_max, n_xy)
            x_coords.append(np.sqrt(2*z)*np.cosh(theta))
            y_coords.append(np.sqrt(2*z)*np.sinh(theta))
            z_coords.append(z*np.ones_like(theta))

    X = np.concatenate(x_coords)
    Y = np.concatenate(y_coords)
    Z = np.concatenate(z_coords)

    # Extend to all 4 quadrants
    X, Y, Z = np.concatenate([X, -X, Y, Y]), np.concatenate([Y, Y, X, -X]), np.concatenate([Z, Z, -Z, -Z])

    # Normal direction at each point:
    normal = np.stack([X, -Y, -np.ones_like(Z)], axis=-1)  # Normal: A x_0 + B -> [x_0, -y_0, -1]
    normal /= np.linalg.norm(normal, axis=-1)[..., np.newaxis]

    # Generate polyhedrons for OpenSCAD
    #   Each point is a list of 3 coordinates
    outer = np.stack([X, Y, Z], axis=-1) + normal * (width_mm / 10) / 2
    inner = outer - normal * (width_mm / 10)
    points_cm = np.concatenate([outer, inner], axis=0)
    # Correct the ranges due to rounding errors, also normal vector can have contributions in each direction
    points_cm[:, 0] = np.maximum(points_cm[:, 0], range_x_cm[0])
    points_cm[:, 0] = np.minimum(points_cm[:, 0], range_x_cm[1])
    points_cm[:, 1] = np.maximum(points_cm[:, 1], range_y_cm[0])
    points_cm[:, 1] = np.minimum(points_cm[:, 1], range_y_cm[1])
    points_cm[:, 2] = np.maximum(points_cm[:, 2], range_z_cm[0])
    points_cm[:, 2] = np.minimum(points_cm[:, 2], range_z_cm[1])
    points_mm = 10*points_cm
    plot(outer, inner)

    #   Each face is a collection of 6 points:
    #       bot_x_first, bot_x_second, bot_y_second, bot_y_first
    #       top_x_first, top_x_second, top_y_second, top_y_first
    idcs = np.arange(n_xy)
    faces_first_layer = np.stack([idcs[1:], idcs[:-1], idcs[:-1] + n_xy, idcs[1:] + n_xy], axis=-1)
    faces_outer_onequadrant = np.concatenate([faces_first_layer + n_xy * i for i in range(n_z - 1)], axis=0)
    faces_outer = np.concatenate([faces_outer_onequadrant + n_xy * n_z * i for i in range(4)], axis=0)
    faces_inner = faces_outer + n_xy * n_z * 4

    faces_join_bottom_one_quadrant = np.stack([idcs[1:], idcs[:-1], idcs[:-1] + 4*n_xy * n_z, idcs[1:] + 4*n_xy * n_z], axis=-1)
    faces_join_top_one_quadrant = faces_join_bottom_one_quadrant + n_xy * (n_z - 1)
    faces_join_top = np.concatenate([faces_join_top_one_quadrant + n_xy * n_z * i for i in range(4)], axis=0)

    idcs2 = np.arange(0, n_xy * n_z, n_xy)
    faces_left_onequadrant = np.stack([idcs2[1:], idcs2[:-1], idcs2[:-1] + 4* n_xy * n_z, idcs2[1:] + 4* n_xy * n_z], axis=-1)
    faces_left = np.concatenate([faces_left_onequadrant + n_xy * n_z * i for i in range(4)], axis=0)
    faces_right = faces_left + n_xy - 1

    faces_edges = np.concatenate([faces_join_top, faces_left, faces_right], axis=0)
    faces = np.concatenate([faces_outer, faces_inner, faces_edges], axis=0)

    # Check points are within the range
    assert np.all((range_x_cm[0] <= points_cm[:, 0]) & (points_cm[:, 0] <= range_x_cm[1])), 'Out of bounds'
    assert np.all((range_y_cm[0] <= points_cm[:, 1]) & (points_cm[:, 1] <= range_y_cm[1])), 'Out of bounds'
    assert np.all((range_z_cm[0] <= points_cm[:, 2]) & (points_cm[:, 2] <= range_z_cm[1])), 'Out of bounds'

    np.set_printoptions(threshold=sys.maxsize)
    os.makedirs(f"scad_{width_mm:.01f}mm", exist_ok=True)
    with open(f"scad_{width_mm:.01f}mm/hyperbolycparaboloid.scad", "w") as f:
        f.write(f"points = {np.array2string(points_mm, separator=', ')};\n\n")
        f.write(f"faces = {np.array2string(faces, separator=', ')};\n\n")
        f.write("polyhedron(points, faces, convexity=1);\n\n")




if __name__ == "__main__":
    generate()