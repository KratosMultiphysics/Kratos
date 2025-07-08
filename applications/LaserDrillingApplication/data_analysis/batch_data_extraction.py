"""
Extract data from a set of files
"""

import argparse
from pathlib import Path
from data_manipulation_functions import read_single_bore, clean_data, remove_surface_and_outliers, subsample_data, calculate_slices

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
    This script takes as its argument the path of a folder containing multiple bore measurements, extracts some information from each dataset and
    exports it to a json
    """
    # ==== Configurable Parameters ====
    sample_x_min = 0
    sample_x_max = 50
    sample_y_min = 0
    sample_y_max = 50
    sample_limits = (sample_x_min, sample_x_max, sample_y_min, sample_y_max)

    max_points = 5000  # Max number of points to use for entire dataset
    subsample_points = False  # Toggle for random subsampling of entire dataset

    slice_thickness = 1.5 # Thickness of each slice in um
    deep_outlier_quantile_threshold = (
        0.0005  # Lower quantile for discarding deep outliers. Set it to None to not discard any.
    )

    shallow_z_threshold = (
        -1
    )  # Discard points shallower than this z (closer to surface). Set it to None to not discard any.
    random_seed = 42  # Seed for reproducibility

    spline_order = 3  # Spline interpolation order
    n_spline_points = 500  # Number of points in the splines



    # Read folder path as argument
    parser = argparse.ArgumentParser(description="Process multiple bore measurements.")
    parser.add_argument("--path", required=True, help="Path to folder containing the data files.")
    parser.add_argument("--material", required=True, help="Material name to associate with all experiments")

    args = parser.parse_args()
    folderpath = args.path
    material = args.material

    folder = Path(folderpath)

    for file in folder.glob("*.dat"):
        filename = file.stem
        print(filename)
        
        power, pulses, identifier = parse_filename(file)
        
        experiment = {
        "power": power,
        "pulses": pulses,
        "id": identifier,
        "material": material,
        "filename": filename
        }

        # ==== Read and Clean the Data ====
        data_raw = read_single_bore(file)
        data = clean_data(data_raw)
        data, surface_outliers, deep_outliers = remove_surface_and_outliers(data, shallow_z_threshold, deep_outlier_quantile_threshold)

        # ==== Random subsample after filtering ====
        if subsample_points and max_points < len(data):
            data = subsample_data(data, max_points, random_seed)

        # ==== Compute data slices ====
        slices, slice_bounds = calculate_slices(data, slice_thickness)

        print(slice_bounds)




