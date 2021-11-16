from os import path, listdir, getcwd

_list_extension = ".post.lst"
_ascii_extension = ".post.res"
_binary_extension = ".post.bin"

def GenerateGiDListFile(first_dir = None, other_dirs = [], file_name = None):
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
    file_names = [f for f in listdir(directory) if path.isfile(path.join(directory, f))]
    post_files = [directory + f for f in file_names if _IsGiDFile(f)]
    return post_files


def _IsGiDFile(file_name):
    if file_name.endswith(_binary_extension):
        return True
    elif file_name.endswith(_ascii_extension):
        return True
    else:
        return False


def _CheckFileName(file_name):
    if not file_name:
        file_name = path.basename(getcwd())
    if not file_name.endswith(_list_extension):
        file_name += _list_extension
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
