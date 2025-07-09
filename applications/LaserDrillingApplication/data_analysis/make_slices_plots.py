import json
import argparse

if __name__ == "__main__":
    """
    Reads ellipses from a file and plots, for each experiment, the grid of slices and their ellipses.
    """
    # Read arguments
    parser = argparse.ArgumentParser(description="Save a plot for each measurement containing slice views of the data and the best fitting ellipse for each slice.")

    parser.add_argument("--input-ellipses", required=True, help="Path to a file containing the ellipses.")
    parser.add_argument("--input-folder", required=True, help="Path of the folder that contains the experimental data.")
    parser.add_argument("--output-path", required=True, help="Path to the output folder were the plots will be saved.")

    args = parser.parse_args()

    ellipses_path = args.input_ellipses
    data_path = args.input_folder
    output_path = args.output_path

    with open(ellipses_path, "r") as ellipses_file:
        experiments = json.load(ellipses_file)
    
