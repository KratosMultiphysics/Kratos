"""
Extract ellipses from a folder with files of geometry point clouds
"""

import argparse
from pathlib import Path
import json
from data_manipulation_functions import (
    read_single_bore,
    clean_data,
    remove_surface_and_outliers,
    subsample_data,
    calculate_slice_bounds,
    calculate_slices,
    compute_ellipses,
)
from plotting_functions import plot_individual_slices_grid


def parse_filename(file):
    """Extract power, pulses, and ID from filename like '10W20P03.dat'"""
    base = file.stem
    w_index = base.index("W")
    p_index = base.index("P")

    power = int(base[:w_index])
    pulses = int(base[w_index + 1 : p_index])
    identifier = base[p_index + 1 :]

    return power, pulses, identifier


if __name__ == "__main__":
    """
    This script takes as its argument the path of a folder containing multiple bore measurements, extracts the ellipses from each dataset and exports them to a json
    """

    # Read arguments
    parser = argparse.ArgumentParser(description="Process multiple bore measurements to extract ellipses.")

    parser.add_argument("--input-folder", required=True, help="Path to a folder containing the data files.")
    parser.add_argument(
        "--parameters-file", required=True, help="Path of the file that contains the parameters to be used"
    )

    parser.add_argument(
        "--output-file", required=True, help="Path to the output JSON file (note that it will be overwritten)."
    )
    parser.add_argument(
        "--output-failures",
        required=True,
        help="Path to the output JSON file that conatins the list of data files that failed and the reason.",
    )
    parser.add_argument(
        "--output-plots",
        required=False,
        help="Path to the output folder were the plots will be saved. If it is not specfied, plots are not generated.",
    )

    parser.add_argument("--material", required=True, help="Material name to associate with all experiments")

    args = parser.parse_args()

    folderpath = args.input_folder
    parameters_file = args.parameters_file

    output_file = args.output_file
    output_failures = args.output_failures
    output_plots_path = args.output_plots

    material = args.material

    input_data_folder = Path(folderpath)

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
        slice_thickness = params["slice_thickness"]  # Thickness of each slice in um
        deep_outlier_quantile_threshold = params[
            "deep_outlier_quantile_threshold"
        ]  # Lower quantile for discarding deep outliers. Set it to None to not discard any.

        shallow_z_threshold = params[
            "shallow_z_threshold"
        ]  # Discard points shallower than this z (closer to surface). Set it to None to not discard any.
        random_seed = params["random_seed"]  # Seed for reproducibility
        spline_order = params["spline_order"]  # Spline interpolation order
        n_spline_points = params["n_spline_points"]  # Number of points in the splines
    except KeyError:
        print("Missing parameters in the parameters JSON file.")
        exit(-1)

    failed_files = {}  # Dictionary of files that failed for some reason
    number_of_files = sum(1 for _ in input_data_folder.glob("*.dat"))
    print(f"{number_of_files=}")

    experiments = {}
    for i, file in enumerate(input_data_folder.glob("*.dat")):
        filename = file.name
        print(f"Analyzing file {filename}: {i + 1}/{number_of_files}")

        power, pulses, identifier = parse_filename(file)

        # ==== Read and Clean the Data ====
        data_raw = read_single_bore(file)
        data = clean_data(data_raw)

        # ALIGN THE DATA TO THE XY PLANE

        data, surface_outliers, deep_outliers = remove_surface_and_outliers(
            data, shallow_z_threshold, deep_outlier_quantile_threshold
        )

        # ==== Random subsample after filtering ====
        if subsample_points and max_points < len(data):
            data = subsample_data(data, max_points, random_seed)

        # ==== Compute data slices ====
        try:
            slice_bounds = calculate_slice_bounds(data, slice_thickness)
        except IndexError:
            failed_files[filename] = "Empty data list"
            print("Error in file" + filename)
            continue

        # plot_cloud_and_slices(x, y, z, slice_bounds, sample_limits, filename)
        try:
            slices = calculate_slices(data, slice_bounds)
        except ValueError:
            # Emtpy slice, handle gracefully so that the script can
            # move on to the next data file. Note the file that is
            # failing for manual review later.
            failed_files[filename] = "Empty slice"
            print("Error in file" + filename)
            continue

        try:
            ellipses = compute_ellipses(slices, slice_bounds)
        except ValueError:
            failed_files[filename] = "Error calculating ellipses"
            print("Error in file" + filename)
            continue

        experiment = {"material": material, "power": power, "pulses": pulses, "id": identifier, "ellipses": ellipses}

        experiments[filename] = experiment

        # Make the plot and export it
        if output_plots_path is not None:
            plot_path = Path(output_plots_path, file.stem)
            plot_individual_slices_grid(
                slices,
                slice_bounds,
                centroids=None,
                ellipses=ellipses,
                sample_limits=sample_limits,
                plot_centroids=False,
                plot_ellipses_in_slices=True,
                filename=filename,
                save_path=plot_path,
            )

    # Export the JSON with the ellipses
    with open(output_file, "w") as fout:
        json.dump(experiments, fout, indent=2)

    # Export the JSON with the failures
    if len(failed_files) > 0:
        with open(output_failures, "w") as fout:
            json.dump(failed_files, fout, indent=2)

        print("The following files failed:" + str(failed_files))
    else:
        print("No files failed")
