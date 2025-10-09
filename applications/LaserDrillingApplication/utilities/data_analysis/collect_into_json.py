import os
import json
import argparse


def parse_filename(filename):
    """Extract power, pulses, and ID from filename like '10W20P03.dat'"""
    base = os.path.splitext(filename)[0]
    w_index = base.index("W")
    p_index = base.index("P")

    power = int(base[:w_index])
    pulses = int(base[w_index + 1 : p_index])
    identifier = base[p_index + 1 :]

    return power, pulses, identifier


def load_point_cloud(filepath):
    """Load the .dat file, returning a list of ['x', 'y', 'z'] as strings"""
    points = []
    with open(filepath, "r") as f:
        for i, line in enumerate(f):
            parts = line.strip().split()
            if len(parts) == 3:
                points.append(parts)
            else:
                raise ValueError("The point in line {i} couldn't be parsed into three values for x, y and z")
    return points


def process_data_folder(folder_path, material):
    experiments = []

    for filename in os.listdir(folder_path):
        if filename.endswith(".dat"):
            power, pulses, identifier = parse_filename(filename)
            filepath = os.path.join(folder_path, filename)
            points = load_point_cloud(filepath)

            experiment = {
                "power": power,
                "pulses": pulses,
                "id": identifier,
                "material": material,
                "filename": filename,
                "points": points,
            }

            experiments.append(experiment)

    return experiments


def save_json(data, output_path):
    with open(output_path, "w") as f:
        json.dump(data, f, indent=2)


def main():
    parser = argparse.ArgumentParser(description="Merge experimental .dat files into a JSON dataset.")
    parser.add_argument("--folder", required=True, help="Path to folder containing .dat files")
    parser.add_argument("--output", required=True, help="Output JSON file path")
    parser.add_argument("--material", required=True, help="Material name to associate with all experiments")

    args = parser.parse_args()

    data = process_data_folder(args.folder, args.material)
    save_json(data, args.output)
    print(f"Saved {len(data)} experiments to '{args.output}' with material '{args.material}'")


if __name__ == "__main__":
    main()
