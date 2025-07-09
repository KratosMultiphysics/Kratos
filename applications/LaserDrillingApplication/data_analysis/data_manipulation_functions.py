"""
Functions to load the data, clean it and fit it
"""

import numpy as np
from scipy.interpolate import splprep, splev
from skimage.measure import EllipseModel
import csv



# ==== Ellipse Fitting Function ====
def fit_ellipse_least_squares(points_2d):
    """
    Fit an ellipse to all given 2D points using least-squares method.

    Parameters:
    - points_2d: Nx2 array of (x, y) points.

    Returns:
    - Tuple: (cx, cy, a, b, eccentricity, theta) or None if fitting fails.
    """
    if len(points_2d) < 5:
        raise ValueError("Not enough points to fit the ellipse")

    model = EllipseModel()
    if not model.estimate(points_2d):
        raise ValueError("Couldn't fit the ellipse to the data")

    cx, cy, a, b, theta = model.params  # (x, y, a, b, theta)
    if b > a:
        a, b = b, a
        theta += np.pi / 2

    eccentricity = np.sqrt(1 - (b**2 / a**2))
    return (cx, cy, a, b, eccentricity, theta)





def calculate_chord_length(a):
    """
    Calculate the chord length of a curve, that is, the distance from its endpoints,
    given by points in array a of shape (N,3) where the columns represent
    the (x,y,z) coordinates of the points
    """
    ax, ay, az = a[:, 0], a[:, 1], a[:, 2]
    return np.linalg.norm([ax[-1] - ax[0], ay[-1] - ay[0], az[-1] - az[0]])


def calculate_arc_length(a):
    """
    Calculate the arc-length of a curve given by points in array a of shape (N,3)
    where the columns represent the (x,y,z) coordinates of the points
    """
    return np.sum(np.linalg.norm(np.diff(np.array(a), axis=0), axis=1))


def calculate_tortuosity(a):
    """
    Calculate the tortuosity of a curve given by points in array a of shape (N,3)
    where the columns represent the (x,y,z) coordinates of the points
    """
    arc_length = calculate_arc_length(a)
    chord_length = calculate_chord_length(a)
    return arc_length / chord_length


def fit_spline_3d(x, y, z, t_arr, spline_order=3):
    """
    Fits a spline through 3D points defined by x, y, z coordinates.

    Parameters:
        x, y, z (array-like): Coordinates of the points.
        t_arr (array of float): Values of the spline parameter where to evaluate the spline.
        spline_order (int): The order of the spline (default: 3).

    Returns:
        spline_points (np.ndarray): Array of shape (num_points, 3) with the fitted spline.
        tck (tuple): The B-spline representation returned by splprep.

    Raises:
        ValueError: If fewer than 2 points are provided or fitting fails.
    """
    if len(x) < 2:
        raise ValueError(f"Need at least 2 points to fit a spline, got {len(x)}.")

    try:
        k = min(spline_order, len(x) - 1)
        tck, _ = splprep([x, y, z], s=0, k=k)
        spline = splev(t_arr, tck)
        spline_points = np.array(spline).T
        return spline, spline_points, tck
    except Exception as e:
        raise ValueError(f"Failed to fit spline: {str(e)}")


def compute_derivatives_spline(t_arr, tck):
    """
    Given a list of points and the knots of a spline, it computes the first,
    second and third derivatives of the spline at said points.
    """
    dx, dy, dz = splev(t_arr, tck, der=1)
    ddx, ddy, ddz = splev(t_arr, tck, der=2)
    dddx, dddy, dddz = splev(t_arr, tck, der=3)

    # Each row represents the (n-th order) derivative at a point in the spline (column 1 = x, column 2 = y, column 3 = z)
    D1 = np.vstack((dx, dy, dz)).T
    D2 = np.vstack((ddx, ddy, ddz)).T
    D3 = np.vstack((dddx, dddy, dddz)).T

    return D1, D2, D3


def compute_curvature(D1, D2):
    """
    Computes the local curvature of a curve given its first (D1) and second (D2) derivatives
    """

    cross = np.cross(D1, D2)
    numerator = np.linalg.norm(cross, axis=1)
    denominator = np.linalg.norm(D1, axis=1) ** 3
    curvature = numerator / denominator
    curvature = np.where(np.isfinite(curvature), curvature, 0)  # Handle division by zero

    return curvature


def compute_torsion(D1, D2, D3):
    """
    Computes the local torsion of a curve given its first (D1), second (D2) and third (D3) derivatives
    """
    cross = np.cross(D1, D2)
    torsion_numerator = np.einsum("ij,ij->i", cross, D3)
    torsion_denominator = np.linalg.norm(cross, axis=1) ** 2
    torsion = torsion_numerator / torsion_denominator
    torsion = np.where(np.isfinite(torsion), torsion, 0)  # Handle division by zero

    return torsion


def read_single_bore(filepath):
    """
    Reads the data of a single bore (from files like those named 10W5P02.dat)
    """
    data = []
    try:
        with open(filepath, "r") as file:
            reader = csv.reader(file, delimiter=" ")
            for row in reader:
                try:
                    point = [float(coord) for coord in row if coord.strip()]
                    if len(point) == 3 and all(np.isfinite(point)):
                        data.append(point)
                except ValueError:
                    continue
    except FileNotFoundError:
        raise FileNotFoundError(f"Error: File {filepath} not found.")

    data = np.array(data)

    return data


def clean_data(data):
    """
    Removes invalid data
    """
    total_points = len(data)

    # Remove non-numerical or NaN entries
    valid_mask = np.isfinite(data).all(axis=1)
    data_clean = data[valid_mask]
    num_nan = total_points - len(data)
    print(f"Discarded {num_nan} non-numerical/NaN points ({100 * num_nan / total_points:.2f}%)")

    return data_clean


def remove_surface_and_outliers(data, shallow_z_threshold, deep_outlier_quantile_threshold):
    """
    Removes the points corresponding to the surface by removing the points shallower than
    shallow_z_threshold (points closer to the surface than shallow_z_threshold), and removes
    outliers that lie too deep in the sample that correspond to wrongly measured data by
    removing the deep_outlier_quantile_threshold deepest points
    """
    total_points = len(data)

    z_vals = data[:, 2]

    # If shallow_z_threshold is None, do not exclude any points from the surface
    surface_mask = np.zeros_like(z_vals, dtype=bool) if shallow_z_threshold is None else z_vals > shallow_z_threshold

    # Similarly, to not exclude deep outliers
    if deep_outlier_quantile_threshold is None:
        deep_mask = np.zeros_like(z_vals, dtype=bool)
    else:
        z_deep_cutoff = np.quantile(z_vals, deep_outlier_quantile_threshold)
        deep_mask = z_vals < z_deep_cutoff

    keep_mask = ~(surface_mask | deep_mask)

    surface_outliers = data[surface_mask]
    deep_outliers = data[deep_mask]
    data = data[keep_mask]

    num_surface_outliers = len(surface_outliers)
    num_deep_outliers = len(deep_outliers)


    # print(f"{num_surface_outliers=}")
    # print(f"{total_points=}")
    print(f"Discarded {num_surface_outliers} surface outliers({100 * num_surface_outliers / total_points:.2f}%)")
    print(f"Discarded {num_deep_outliers} deep outliers({100 * num_deep_outliers / total_points:.2f}%)")

    return data, surface_outliers, deep_outliers


def subsample_data(data, max_points, random_seed):
    """
    Subsamples the data by randomly keeping max_points and discarding the rest
    """

    np.random.seed(random_seed)
    indices = np.random.choice(len(data), max_points, replace=False)
    data = data[indices]

    return data

def calculate_slice_bounds(data, slice_thickness):
    """
    Calculates the slice bounds
    """
    z = data[:, 2]
    if z.size <= 0:
        raise IndexError("Empty list of z values in the data received by calculate_slice_bounds.")

    # ==== Create Slices Variable ====
    z_min, z_max = z.min(), z.max()
    print(f"z_max={float(z_max)}, z_min={float(z_min)}")


    num_slices = int(np.ceil((z_max - z_min)/slice_thickness))

    # Generate slice boundaries from top (z_max) to bottom (z_min)
    slice_bounds = np.linspace(z_max, z_max - num_slices*slice_thickness, num_slices+1)
    # print(slice_bounds)
    return slice_bounds

def calculate_slices(data, slice_bounds):
    """
    Calculates the data slices
    """
    z = data[:, 2]

    num_slices = slice_bounds.size-1
    # Initialize list to store slices
    slices = []

    # Create slices
    for i in range(num_slices):
        z_start, z_end = slice_bounds[i], slice_bounds[i + 1]
        # print(f"{i=}")
        # print(f"{z_start=}, {z_end=}")
        
        # Mask for points in the current slice
        mask = (z <= z_start) & (z > z_end)
        
        data_in_slice = data[mask]
        # print(data_in_slice)

        if len(data_in_slice) <= 0:
            raise ValueError(f"Empty slice between {z_start=}, {z_end=}")

        slices.append(data_in_slice)

    return slices

def compute_centroids(slices, slice_bounds):
    """
    For each slice, finds its centroid
    """
    num_slices = len(slices)

    centroids = [None] * num_slices
    centroid_zs = [None] * num_slices
    centroid_ids = [None] * num_slices

    for i in range(num_slices):
        print(f"Computing centroid for Slice {i + 1} (z ∈ [{slice_bounds[i]:.2f}, {slice_bounds[i + 1]:.2f}] um)")
        slice_points = slices[i]

        if len(slice_points) == 0:
            raise ValueError(
                f"Slice {i + 1} (z ∈ [{slice_bounds[i]:.2f}, {slice_bounds[i + 1]:.2f}] um) contains no points"
            )

        centroid = np.mean(slice_points, axis=0)

        centroids[i] = centroid
        centroid_ids[i] = i + 1

    # Filter out None values
    valid_centroids = [c for c in centroids if c is not None]
    centroids = np.array(valid_centroids)
    centroid_zs = [z for z in centroid_zs if z is not None]
    centroid_ids = [id for id in centroid_ids if id is not None]

    return centroids, centroid_ids


def compute_ellipses(slices, slice_bounds):
    """
    Find the ellipse that best fits each spline using least squares
    """
    num_slices = len(slices)

    ellipses = [None] * num_slices

    for i in range(num_slices):
        print(f"Fitting ellipse for Slice {i + 1} (z ∈ [{slice_bounds[i]:.2f}, {slice_bounds[i + 1]:.2f}] um)")
        slice_points = slices[i]

        if len(slice_points) == 0:
            raise ValueError("Empty slice")

        points_2d = slice_points[:, :2]  # [x, y]  

        center_x, center_y, a, b, eccentricity, angle = fit_ellipse_least_squares(points_2d)
        center_z = np.mean(slice_points[:, 2]) # Set the ellipse's depth as the "center of mass" in depth of the points in the slice

        ellipse_center = np.array([center_x, center_y, center_z])

        ellipse_dict = {"center": ellipse_center, "a": a, "b": b, "eccentricity": eccentricity, "angle": angle}
        ellipses[i] = ellipse_dict            
            

    return ellipses
