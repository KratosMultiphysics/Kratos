from data_manipulation_functions import read_single_bore, clean_data, remove_surface_and_outliers, subsample_data, calculate_slice_bounds, calculate_slices, compute_centroids, compute_ellipses, fit_spline_3d, compute_derivatives_spline, compute_curvature, compute_torsion

from plotting_functions import plot_3d_geometry, plot_outliers, plot_xy_centers_path, plot_superposed_slices_xy, plot_ellipses_and_axes, plot_contour, plot_ellipse_metrics, plot_individual_slices_grid, plot_individual_slices_separately

import numpy as np
import argparse
from pathlib import Path
import json

if __name__ == "__main__":
    """
    This script takes as its argument a name of a dataset corresponding to one bore and makes multiple plots.
    The default actions leave the data unchanged. Any modifications to the data, such as removing outliers or rotating it must be toggled on by the user explicitely.
    """

    # Read filepath as argument
    parser = argparse.ArgumentParser(description="Make plots and extract metrics of the measurements of a bore.")


    parser.add_argument("--filepath", required=True, help="Path to file containing the data file.")

    parser.add_argument("--parameters-file", required=True, help="Path of the file that contains the parameters to be used")
    
    # The sample surface won't exist if it has been removed. Therefore, these arguments can't be both true simultaneously
    arg_group = parser.add_mutually_exclusive_group(required=False)
    arg_group.add_argument("--align-data", action="store_true", help="Toggle to rotate and translate the data points so that the sample's face is coincident with the XY plane. Incompatible with --remove-surface-and-outliers")

    arg_group.add_argument("--remove-surface-and-outliers", action="store_true", help="Toggle to remove points above the surface and outliers based on depth. Optional (by default, they are NOT removed). Incompatible with --align-data")

    args = parser.parse_args()

    filepath = args.filepath
    parameters_file = args.parameters_file

    align_data_toggle = args.align_data
    remove_surface_and_outliers_toggle = args.remove_surface_and_outliers
    

    print(f"{remove_surface_and_outliers_toggle=}")
    filename = Path(filepath).stem
    print(f"Plots for {filename}")


        # Read parameters
    with open(parameters_file, "r") as params_file:
        params = json.load(params_file)

    # ==== Configure Parameters ====
    try:
        sample_x_min = params["sample_x_min"]
        sample_x_max = params["sample_x_max"]
        sample_y_min = params["sample_y_min"]
        sample_y_max = params["sample_y_max"]
        sample_limits = (sample_x_min, sample_x_max, sample_y_min, sample_y_max)

        max_points = params["max_points"]  # Max number of points to use for entire dataset
        subsample_points = params["subsample_points"]  # Toggle for random subsampling of entire dataset
        slice_thickness = params["slice_thickness"] # Thickness of each slice in um
        deep_outlier_quantile_threshold = params["deep_outlier_quantile_threshold"]  # Lower quantile for discarding deep outliers. Set it to None to not discard any.
        
        shallow_z_threshold = params["shallow_z_threshold"] # Discard points shallower than this z (closer to surface). Set it to None to not discard any.
        random_seed = params["random_seed"]  # Seed for reproducibility
        spline_order = params["spline_order"]  # Spline interpolation order
        n_spline_points = params["n_spline_points"]  # Number of points in the splines
    except KeyError:
        print("Missing parameters in the parameters JSON file.")
        exit(-1)

    # Plot toggles
    plot_3d_geometry_toggle = True  # Plot of the point cloud with the fitted spline(s) through the center(s)
    plot_planes = False
    plot_outliers_toggle = False
    plot_outlier_types = "deep"  # Options: "all", "surface", "deep"
    plot_xy_center_approximations = True
    plot_xy_centroids = False
    plot_xy_ellipse_centers = True
    plot_z_histogram = False  # Optional histogram to guide shallow_z_threshold
    plot_superposed_slices_xy_toggle = True  # Toggle for XY view of all slices superposed
    plot_centroids_in_superposed_slices_xy = False  # Toggle for centroids in XY view of superposed slices
    plot_ellipses_in_superposed_slices_xy = True  # Toggle for centroids in XY view of superposed slices
    plot_ellipse_centers_in_superposed_slices_xy = True  # Toggle for ellipses' centers in XY view of superposed slices
    plot_ellipses_and_axes_toggle = True  # Toggle for plot the of all ellipses
    ellipses_shift_to_common_center = False  # Toggle to shift all ellipses so that they have the same center
    plot_contour_toggle = False  # Toggle for contour plot
    plot_ellipses_in_slices = True  # Toggle for plotting ellipses in slice views
    plot_ellipse_metrics_toggle = True  # Toggle for eccentricity and theta vs. depth plot
    plot_individual_slices_grid_toggle = True
    plot_individual_slices_grid_centroids = False
    plot_slice_views_one_by_one = False  # Toggle for XY view of each slices

    # ==== Read and Clean the Data ====
    data_raw = read_single_bore(filepath)
    data = clean_data(data_raw)

    if remove_surface_and_outliers_toggle:
        data, surface_outliers, deep_outliers = remove_surface_and_outliers(data, shallow_z_threshold, deep_outlier_quantile_threshold)

    # ==== Random subsample after filtering ====
    if subsample_points and max_points < len(data):
        data = subsample_data(data, max_points, random_seed)

    x, y, z = data[:, 0], data[:, 1], data[:, 2]

    # ==== Compute data slices ====
    try:
        slice_bounds = calculate_slice_bounds(data, slice_thickness)
    except IndexError:
        print("Empty data list")
        exit(-1)

    plot_3d_geometry(x, y, z, sample_limits, filename, slice_bounds=slice_bounds if plot_planes else None)

    try:
        slices = calculate_slices(data, slice_bounds)
    except ValueError:
        print("There are empty slices")
        exit(-1)


    # ==== Compute centroids for slices ====
    centroids, centroid_ids = compute_centroids(slices, slice_bounds)

    centroids_x, centroids_y, centroids_z = centroids[:, 0], centroids[:, 1], centroids[:, 2]

    # ==== Fit ellipses to slices ====
    try:
        ellipses = compute_ellipses(slices, slice_bounds)
    except ValueError:
        print("Errors computing the ellipses")
        exit(-1)


    ellipse_centers = []
    for ellipse in ellipses:
        ellipse_centers.append(ellipse["center"])

    ellipse_centers = np.array(ellipse_centers)

    ellipses_x, ellipses_y, ellipses_z = ellipse_centers[:, 0], ellipse_centers[:, 1], ellipse_centers[:, 2]

    # ==== Fit splines ====
    t_arr = np.linspace(0, 1, n_spline_points)

    centroid_spline, centroid_spline_np, centroid_tck = fit_spline_3d(
        centroids_x, centroids_y, centroids_z, t_arr, spline_order=spline_order
    )

    ellipse_spline, ellipse_spline_np, ellipse_tck = fit_spline_3d(
        ellipses_x, ellipses_y, ellipses_z, t_arr, spline_order=spline_order
    )

    # ==== Compute Curvature and Torsion (for Centroid Spline) ====
    D1_centroid, D2_centroid, D3_centroid = compute_derivatives_spline(t_arr, centroid_tck)
    curvature_centroid = compute_curvature(D1_centroid, D2_centroid)
    torsion_centroid = compute_torsion(D1_centroid, D2_centroid, D3_centroid)

    avg_curvature = np.mean(curvature_centroid)
    avg_torsion = np.mean(torsion_centroid)
    print(f"Average Curvature: {avg_curvature:.6f} um^-1")
    print(f"Average Torsion: {avg_torsion:.6f} um^-1")

    # ==== Plotting ====

    if plot_3d_geometry_toggle:
        plot_3d_geometry(
            x,
            y,
            z,
            sample_limits,
            filename,
            centroid_spline,
            centroids_x,
            centroids_y,
            centroids_z,
            ellipse_spline,
            ellipses_x,
            ellipses_y,
            ellipses_z,
            slice_bounds=slice_bounds if plot_planes else None            
        )

    if plot_outliers_toggle:
        plot_outliers(x, y, z, surface_outliers, deep_outliers, plot_outlier_types, filename)

    if plot_xy_center_approximations:
        plot_xy_centers_path(
            plot_xy_centroids,
            plot_xy_ellipse_centers,
            centroids_x,
            centroids_y,
            centroids_z,
            centroid_ids,
            ellipse_centers,
            sample_limits,
            filename
        )

    if plot_superposed_slices_xy_toggle:
        plot_superposed_slices_xy(
            slices,
            centroids,
            centroid_ids,
            ellipses,
            sample_limits,
            plot_centroids_in_superposed_slices_xy,
            plot_ellipses_in_all_slices_xy=True,
            plot_ellipse_centers=True,
            filename=filename
        )

    if plot_ellipses_and_axes_toggle:
        plot_ellipses_and_axes(
            ellipses,
            sample_limits,
            slice_bounds,
            filename,
            shift_to_common_center=ellipses_shift_to_common_center
        )

    if plot_ellipses_and_axes_toggle:
        plot_ellipses_and_axes(
            ellipses,
            sample_limits,
            slice_bounds,
            filename,
            shift_to_common_center=not ellipses_shift_to_common_center
        )

    if plot_contour_toggle:
        contour_levels = len(slices)

        plot_contour(x, y, z, contour_levels, sample_limits, filename)

    if plot_ellipse_metrics_toggle:
        plot_ellipse_metrics(ellipses, filename)

    if plot_individual_slices_grid_toggle:
        plot_individual_slices_grid(slices, slice_bounds, centroids, ellipses, sample_limits, plot_individual_slices_grid_centroids, plot_ellipses_in_slices, filename)

    if plot_slice_views_one_by_one:
        plot_individual_slices_separately(
            slices,
            slice_bounds,
            centroids,
            centroid_ids,
            ellipses,
            sample_limits,
            plot_ellipses_in_slices,
            filename
        )
