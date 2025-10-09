"""
Functions to load the data, clean it and fit it
"""

import numpy as np

from scipy.interpolate import splprep, splev
from scipy.spatial.transform import Rotation
from scipy.spatial.transform import RigidTransform as Tf

from skimage.measure import EllipseModel, ransac

from sklearn.linear_model import RANSACRegressor, LinearRegression

import csv

from pathlib import Path


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


def fit_ellipse_ransac(points_2d, min_samples, residual_threshold, max_trials, rng_seed):
    """
    Fit an ellipse to all given 2D points using RANSAC method.


    See: https://scikit-image.org/docs/stable/api/skimage.measure.html#skimage.measure.EllipseModel

    See the examples in https://scikit-image.org/docs/stable/api/skimage.measure.html#skimage.measure.ransac

    Parameters:
    - points_2d: Nx2 array of (x, y) points.


    Returns:
    - Tuple: (cx, cy, a, b, eccentricity, theta) or None if fitting fails.

    - inlier_mask (np.ndarray): Boolean mask of inliers classified as True.
    """

    try:
        ransac_model, inlier_mask = ransac(
            points_2d, EllipseModel, min_samples=min_samples, residual_threshold=residual_threshold, max_trials=max_trials, rng=rng_seed
        )
    except ValueError as e:
        print(e)
        raise ValueError from e

    cx, cy, a, b, theta = ransac_model.params
    if b > a:
        a, b = b, a
        theta += np.pi / 2

    eccentricity = np.sqrt(1 - (b**2 / a**2))

    return (cx, cy, a, b, eccentricity, theta), inlier_mask


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


def clean_and_export_csv(input_path, output_path=None, fmt=["%.3f", "%.3f", "%.6f"], delimiter=","):
    """
    Reads the data file from a bore, removes the invalid data and exports the clean data to a csv.
    Useful to visualize the files in e.g. Paraview
    """

    input_path = Path(input_path)
    print("Cleaning and exporting as CSV" + input_path.name)

    data_clean = clean_data(read_single_bore(input_path))

    if output_path is None:
        output_path = input_path.with_stem(input_path.stem + "_clean").with_suffix(".csv")

    try:
        np.savetxt(output_path, data_clean, fmt, delimiter)
    except:
        raise
    else:
        print("Exported as " + output_path)


def remove_surface(data, shallow_z_threshold):
    """
    Removes the points corresponding to the surface by removing the points shallower than
    shallow_z_threshold (points closer to the surface than shallow_z_threshold)

    Parameters
    ----------
    data (ndarray): Nx3 array with data points
    shallow_z_threshold (float): points with z larger than shallow_z_threshold are considred to be part of the surface and are removed

    Returns
    -------
    data (ndarray): the data without the surface points
    surface (ndarray): points corresponding to the surface
    """

    total_points = len(data)

    z_vals = data[:, 2]

    surface_mask = z_vals > shallow_z_threshold

    keep_mask = ~surface_mask

    surface = data[surface_mask]
    data = data[keep_mask]

    n_surface = len(surface)
    print(f"Discarded {n_surface} surface outliers({100 * n_surface / total_points:.2f}%)")

    return data, surface


def remove_deep_outliers(data, quantile_threshold):
    """
    Removes outliers that lie too deep in the sample that correspond to wrongly measured data by removing the deep_outlier_quantile_threshold deepest points

    Parameters
    ----------
    data (ndarray): Nx3 array with data points
    quantile_threshold (float): points whose z is in the quantile determined by quantile_threshold are removed

    Returns
    -------
    data (ndarray): the data without the deep points
    deep_outliers (ndarray): points considered deep outliers
    """
    total_points = len(data)

    z_vals = data[:, 2]

    z_deep_cutoff = np.quantile(z_vals, quantile_threshold)
    deep_mask = z_vals < z_deep_cutoff

    keep_mask = ~deep_mask

    deep_outliers = data[deep_mask]
    data = data[keep_mask]

    num_deep_outliers = len(deep_outliers)

    print(f"Discarded {num_deep_outliers} deep outliers({100 * num_deep_outliers / total_points:.2f}%)")

    return data, deep_outliers



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

    num_slices = int(np.ceil((z_max - z_min) / slice_thickness))

    # Generate slice boundaries from top (z_max) to bottom (z_min)
    slice_bounds = np.linspace(z_max, z_max - num_slices * slice_thickness, num_slices + 1)
    # print(slice_bounds)
    return slice_bounds


def calculate_slices(data, slice_bounds):
    """
    Calculates the data slices
    """
    z = data[:, 2]

    num_slices = slice_bounds.size - 1
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
            # raise ValueError(f"Empty slice between {z_start=}, {z_end=}")
            print(f"Empty slice between {z_start=}, {z_end=}, appending the empty slice to the list anyways")
            
        slices.append(data_in_slice)

    return slices


def compute_centroids(slices, slice_bounds):
    """
    For each slice, finds its centroid
    """
    num_slices = len(slices)

    centroids = [None] * num_slices
    centroid_ids = [None] * num_slices

    for i in range(num_slices):
        print(f"Computing centroid for Slice {i + 1} (z ∈ [{slice_bounds[i]:.2f}, {slice_bounds[i + 1]:.2f}] um)")
        slice_points = slices[i]

        if len(slice_points) == 0:
            # raise ValueError(
            #     f"Slice {i + 1} (z ∈ [{slice_bounds[i]:.2f}, {slice_bounds[i + 1]:.2f}] um) contains no points"
            # )
            print(f"Slice {i + 1} (z ∈ [{slice_bounds[i]:.2f}, {slice_bounds[i + 1]:.2f}] um) contains no points")
            # centroids[i] = None # Not needed since centroids starts full of None
            continue

        centroid = np.mean(slice_points, axis=0)

        centroids[i] = centroid
        centroid_ids[i] = i + 1

    # # Filter out None values
    # valid_centroids = [c for c in centroids if c is not None]
    # centroids_np = np.array(valid_centroids)
    
    # valid_centroid_ids = [id for id in centroid_ids if id is not None]
    # centroid_ids_np = np.array(valid_centroid_ids)
    
    return centroids, centroid_ids


def compute_ellipses(slices, slice_bounds=None, method="least_squares", ransac_params=None, return_inlier_masks=False):
    """
    Find the ellipse that best fits each slice.
    Possible methods for fitting the ellipse are "least_squares" and "ransac".

    Parameters:
    - slices (list): list with all the slices
    - slice_bounds (list): list with the bounds that define the slices. Optional, only used for printing the status.
    - method (string): method to use for the ellipse fitting
    - ransac_params (dict): parameters for the RANSAC method
    - return_inlier_masks (bool): toggle to return the mask of inliers "belonging" in the ellipse

    Returns:
    - ellipses (list): list of the ellipses corresponding to each slice
    - inliers_per_slice (list): if return_inlier_masks is True, return a list where each element is a list of the indices
        of the points that are inliers of the best fit ellipse of the corresponding slice
    """
    print(f"Fitting ellipses using {method}, from the available methods least_squares and ransac.")
    if method == "ransac":
        if ransac_params is None:
            raise ValueError("Parameters for the RANSAC fit need to be specified")
        else:
            try:
                min_samples = ransac_params["min_samples"]
                residual_threshold = ransac_params["residual_threshold"]
                max_trials = ransac_params["max_trials"]
                rng_seed = ransac_params["rng_seed"]
            except KeyError as e:
                raise KeyError("Incorrect parameters for the RANSAC fitting specified") from e

    num_slices = len(slices)

    ellipses = [None] * num_slices
    inliers_per_slice = [None] * num_slices

    for i in range(num_slices):
        if slice_bounds is not None:
            print(f"Fitting ellipse for Slice {i + 1} (z ∈ [{slice_bounds[i]:.2f}, {slice_bounds[i + 1]:.2f}] um)")

        slice_points = slices[i]

        # if len(slice_points) == 0:
        #     raise ValueError("Empty slice")

        # if len(slice_points) < ransac_params["min_samples"]:
        #     print()
        #     ellipses[i] = None
        #     continue

        points_2d = slice_points[:, :2]  # [x, y]

        if method == "least_squares":
            try:
                center_x, center_y, a, b, eccentricity, angle = fit_ellipse_least_squares(points_2d)
            except ValueError as e:
                print("Error fitting ellipse")
                if slice_bounds is not None:
                    print(f"to slice {i + 1} (z ∈ [{slice_bounds[i]:.2f}, {slice_bounds[i + 1]:.2f}] um)")
                print(e)
                ellipses[i] = None
                continue

        elif method == "ransac":
            try:
                (center_x, center_y, a, b, eccentricity, angle), inlier_mask = fit_ellipse_ransac(
                points_2d, min_samples, residual_threshold, max_trials, rng_seed
            )
                inliers_per_slice[i] = inlier_mask
            except ValueError as e:
                print("Error fitting ellipse")
                if slice_bounds is not None:
                    print(f"to slice {i + 1} (z ∈ [{slice_bounds[i]:.2f}, {slice_bounds[i + 1]:.2f}] um)")
                print(e)
                ellipses[i] = None
                continue

        else:
            raise ValueError("The method selected is not available")

        # Set the ellipse's depth as the "center of mass" in depth of the points in the slice
        center_z = np.mean(slice_points[:, 2])

        ellipse_center = [center_x, center_y, center_z]

        ellipse_dict = {"center": ellipse_center, "a": a, "b": b, "eccentricity": eccentricity, "angle": angle}
        ellipses[i] = ellipse_dict

    if method == "least_squares":
        return ellipses
    elif method == "ransac":
        if return_inlier_masks:
            return ellipses, inliers_per_slice
        else:
            return ellipses


def find_sample_face(data, ransac_params=None):
    """
    Finds the face of the sample, that is, the original flat surface before perforation.
    It assumes the data is valid.

    Parameters:
    - data (np.ndarray): Nx3 matrix with the data points
    - ransac_params (dict, default=None): parameters to use for the ransac fitting

    Returns:
    - a, b, c (floats): coefficients of the plane face: z = a*x + b*y + c
    """

    XY = data[:, [0,1]]
    z = data[:, 2]

    # Fit RANSAC with built-in LinearRegression
    if ransac_params is None:
        raise ValueError("Parameters for the RANSAC fit need to be specified")
    else:
        try:
            min_samples = ransac_params["min_samples"]
            residual_threshold = ransac_params["residual_threshold"]
            max_trials = ransac_params["max_trials"]
            rng_seed = ransac_params["rng_seed"]
        except KeyError:
            raise KeyError("Incorrect parameters for the RANSAC fitting specified")
    
    ransac = RANSACRegressor(LinearRegression(), residual_threshold=residual_threshold, min_samples=min_samples, max_trials=max_trials, random_state=rng_seed)
    ransac.fit(XY, z)

    # Coefficients of the plane from LinearRegression: z = a*x + b*y + c.    
    a, b = ransac.estimator_.coef_    
    c = ransac.estimator_.intercept_
    
    return a, b, c

def unit_vector(vector, epsilon=1e-10):
    """
    Returns the normalized vector of the input vector.

    Parameters:
        - vector (ndarray): (N,) vector to be normalized
    
    Returns:
        - (ndarray): (N,) normalized vector
    """
    norm = np.linalg.norm(vector)
    
    if norm < epsilon:
        raise ValueError("Cannot normalize a vector close to zero")
        
    return vector/norm


def angle_between(v1, v2):
    """ 
    Returns the angle in radians between vectors v1 and v2.

    Parameters:
        - v1 (ndarray): (N,) vector.
        - v2 (ndarray): (N,) vector.
    
    Returns:
        - angle (float): angle between v1 and v2 in radians.
    """

    v1_u = unit_vector(v1)
    v2_u = unit_vector(v2)
    return np.arccos(np.clip(np.dot(v1_u, v2_u), -1.0, 1.0))


def rotate_data(data, vector_initial, vector_final):
    """
    Rotates points in data by the rotation that takes initial_vector to final_vector.
    Assumes valid data.
    
    Parameters:
    - data (np.ndarray): Nx3 matrix with the data points.
    - vector_initial (np.ndarray): (3,) array with the initial vector.
    - vector_final (np.ndarray): (3,) array with the final vector.

    Returns:
    - data_rotated (np.ndarray): Nx3 matrix with the rotated data points.
    """

    rot_dir = unit_vector(np.cross(vector_initial, vector_final))

    rot_angle = angle_between(vector_initial, vector_final)

    rot_vector = rot_dir*rot_angle

    r = Rotation.from_rotvec(rot_vector)

    rot_mat = r.as_matrix()

    data_rotated = (rot_mat @ data.T).T

    return data_rotated

def align_data(data, ransac_params=None):
    """
    Aligns the data to the XY plane. It finds the data surface face and rotates and translates the data so that this face
    coincides with the XY plane.

    Parameters:
    - data (np.ndarray): Nx3 matrix with the data points
    - ransac_params (dict, default=None): parameters to use for the ransac fitting

    Returns:
    - data_aligned (np.ndarray): Nx3 matrix with the aligned data points
    """

    # Find the data's plane face
    a, b, c = find_sample_face(data, ransac_params)

    # Rotate the data so it is parallel to the XY plane
    z_dir = np.asarray([0,0,1])

    # Coefficients of the plane from find_sample_face: z = a*x + b*y + c. 
    # Implicit plane equation: a*x + b*y - z + d = 0
    surface_normal_vect = -np.asarray([a, b, -1])

    data_rotated = rotate_data(data, surface_normal_vect, z_dir)

    # Translate the data to the XY plane
    translation_to_XY_vector = -c*z_dir # c is the Z axis intercept of the data's face

    translation_to_XY = Tf.from_translation(translation_to_XY_vector)
    
    data_on_XY = translation_to_XY.apply(data_rotated)

    return data_on_XY
    


def center_data_XY(data, center):
    """
    Centers the data point cloud on the XY plane by translating all points so that their center of mass becomes 'center'

    Parameters
    ----------
    data (ndarray): Nx3 array of points
    center (array): (3,) array with the coordinates of the new center

    Returns
    -------
    data_centered (ndarray): Nx3 array of points whose center of mass on the XY plane is 'center'
    """
    
    # Translate the data on the XY plane to the origin
    com_XY = np.sum(data, axis=0)/data.shape[0]
    com_XY[2] = 0

    translation_to_origin = Tf.from_translation(center-com_XY)

    data_centered = translation_to_origin.apply(data)

    return data_centered