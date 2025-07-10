from data_manipulation_functions import read_single_bore, clean_data
from pathlib import Path
import argparse

if __name__ == "__main__":
    """
    This script finds the size of the sample (xmin, xmax ymin, ymax)
    """

    # Read arguments
    parser = argparse.ArgumentParser(description="Finds the size of the sample")

    parser.add_argument("--input-folder", required=True, help="Path to a folder containing the data files.")

    args = parser.parse_args()
    
    input_data_folder = Path(args.input_folder)

    min_x_dict = {}
    max_x_dict = {}
    min_y_dict = {}
    max_y_dict = {}

    for i, file in enumerate(input_data_folder.glob("*.dat")):
        print(f"{file.stem}  ---  {i}")
        data = clean_data(read_single_bore(file))

        min_x_dict[file.stem] = data[:, 0].min()
        max_x_dict[file.stem] = data[:, 0].max() 

        min_y_dict[file.stem] = data[:, 1].min()
        max_y_dict[file.stem] = data[:, 1].max()

    global_min_x_key = min(min_x_dict, key=min_x_dict.get)
    global_max_x_key = max(max_x_dict, key=max_x_dict.get)

    global_min_y_key = min(min_y_dict, key=min_y_dict.get)
    global_max_y_key = max(max_y_dict, key=max_y_dict.get)

    print("Global results:")
    print(f"min_x={min_x_dict[global_min_x_key]} at {global_min_x_key}")
    print(f"max_x={max_x_dict[global_max_x_key]} at {global_max_x_key}")
    
    print(f"min_y={min_y_dict[global_min_y_key]} at {global_min_y_key}")
    print(f"max_y={max_y_dict[global_max_y_key]} at {global_max_y_key}")



