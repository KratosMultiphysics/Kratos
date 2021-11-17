from pathlib import Path


def GenerateGiDListFile(first_dir = "", other_dirs = [], file_name = ""):
    first_dir = _CheckDirectory(first_dir)
    files = _GetPostFilesList(first_dir)
    for dir in other_dirs:
        dir = _CheckDirectory(dir)
        files += _GetPostFilesList(dir)
    files.sort()
    file_name = _CheckFileName(file_name)
    _WriteListFile(files, file_name)


def _CheckDirectory(directory):
    if not directory:
        directory = "./"
    if not directory.endswith("/"):
        directory += "/"
    return directory


def _GetPostFilesList(directory):
    files = list(Path(directory).glob("*.post.bin"))
    files.extend(list(Path(directory).glob("*.post.res")))
    return files


def _CheckFileName(file_name):
    list_extension = ".post.lst"
    if not file_name:
        file_name = Path.cwd().name
    if not file_name.endswith(list_extension):
        file_name += list_extension
    return file_name


def _WriteListFile(file_list, file_name):
    file = open(file_name, "w")
    file.write("Merge\n")
    for output_file in file_list:
        file.write("{}\n".format(output_file))
    file.close()


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="Generate a GiD .post.lst file to merge the results of several simulations.")
    parser.add_argument("first_location",  nargs="?", type=str, help="Where the post files are located. By default it is the current directory")
    parser.add_argument("other_locations", nargs="*", type=str, help="Add more post files. By default it is empty")
    parser.add_argument("-n", "--name", type=str, help="Name of the output file. By default it is the current folder name")
    args = parser.parse_args()
    GenerateGiDListFile(args.first_location, args.other_locations, args.name)
