def get_data_points_from_file(file_path):
    result = []
    with open(file_path, "r") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue

            _, x, y = line.split()  # ignore the row number
            result.append((float(x), float(y)))

    return result