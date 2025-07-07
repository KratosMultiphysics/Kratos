import os

from bending_measurement_ellipses_simple_encapsulated_plots import read_single_bore, clean_data, remove_surface_and_outliers


def parse_filename(filename):
    """Extract power, pulses, and ID from filename like '10W20P03.dat'"""
    base = os.path.splitext(filename)[0]
    w_index = base.index("W")
    p_index = base.index("P")

    power = int(base[:w_index])
    pulses = int(base[w_index + 1 : p_index])
    identifier = base[p_index + 1 :]

    return power, pulses, identifier


def make_list_of_deepest_point(path):
    """
    Iterates through all files (representing a bore) in a folder and, for each bore,
    it finds the deepest point
    """
    pass

if __name__ == "__main__":
    folder_path ="jkasfhjdsjpafjb"

    for filename in os.listdir(folder_path):
        if filename.endswith(".dat"):
            power, pulses, identifier = parse_filename(filename)
            filepath = os.path.join(folder_path, filename)

            # ==== Read and Clean the Data ====
            data_raw = read_single_bore(filepath)
            data = clean_data(data_raw)
            data, surface_outliers, deep_outliers = remove_surface_and_outliers(
                data, shallow_z_threshold=None, deep_outlier_quantile_threshold
            )