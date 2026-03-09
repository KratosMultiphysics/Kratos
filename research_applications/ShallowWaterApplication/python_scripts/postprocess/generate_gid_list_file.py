from pathlib import Path


def GenerateGiDListFile(directories = [], file_name = ""):
    files = []
    for dir in directories:
        files += _GetPostFilesList(Path(dir))
    files.sort()
    file_name = _ValidateFileName(file_name)
    _WriteListFile(files, file_name)


def _GetPostFilesList(directory):
    files = list(directory.glob("*.post.bin"))
    files.extend(list(directory.glob("*.post.res")))
    return files


def _ValidateFileName(file_name):
    if not file_name:
        file_name = Path.cwd().name
    return Path(file_name).with_suffix("").with_suffix(".post.lst")


def _WriteListFile(file_list, file_name):
    with open(file_name, "w") as file:
        file.write("Merge\n")
        for output_file in file_list:
            file.write("{}\n".format(output_file))


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="Generate a GiD .post.lst file to merge the results from several simulations.")
    parser.add_argument("directory", nargs="*", type=str, default=["."], help="where the post files are located. By default, the first location is the current directory")
    parser.add_argument("-n", "--name", type=str, help="name of the output file. By default it is the current folder name")
    args = parser.parse_args()
    GenerateGiDListFile(args.directory, args.name)
