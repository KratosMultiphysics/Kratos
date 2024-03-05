import os
import fnmatch
from unittest import TestCase
import filecmp
import run_multiple_stages as rms

def detect_subfolders(root_folder):
    """
    Detect every folder inside the specified root folder.

    Parameters:
    - root_folder: The path to the root folder to start the search.

    Returns:
    - A list of subfolder names.
    """
    subfolders = [f.path for f in os.scandir(root_folder) if f.is_dir()]

    return subfolders

def list_files_with_structure(folder, target_structure):
    """
    List files with a specific structure within a given folder.

    Parameters:
    - folder: The path to the folder to search for files.
    - target_structure: The specific structure to match (e.g., file extension).

    Returns:
    - A list of matching files in the specified folder.
    """
    matching_files = [file for file in os.listdir(folder) if fnmatch.fnmatch(file, target_structure)]
    return matching_files

def stages_in_test():
    """
    Detect the number of stages in the current test. i.e. current folder

    Returns:
    - The number of stages in the current test.
    """
    target_file_structure = "ProjectParameters_*"
    matching_files_list = list_files_with_structure(os.getcwd(), target_file_structure)
    if matching_files_list:
        for file in matching_files_list:
            print(f"  File: {file}")
    else:
        raise RuntimeError("No test files found in the current folder.")
    return len(matching_files_list)

def run_test_in_folder():
    """
    Run the test in the current folder.
    """
    print(f"Running test in folder: {os.getcwd()}")
    no_stages = stages_in_test()
    # Run the test
    rms.run_stages(os.getcwd(), no_stages)

    groundtruth_files = list_files_with_structure(os.getcwd(), "*.post.orig.res")
    output_files = list_files_with_structure(os.getcwd(), "*.post.res")

    groundtruth_files.sort()
    output_files.sort()

    if len(groundtruth_files) != len(output_files):
        raise RuntimeError("The number of output files does not match the number of groundtruth files.")

    for truth, output in zip(groundtruth_files, output_files):
        assert compare_two_files(truth, output), f"Output file {output} does not match groundtruth file {truth}."

def compare_two_files(file1, file2):
    """
    Compare two files.

    Parameters:
    - file1: The path to the first file.
    - file2: The path to the second file.

    Returns:
    - A boolean indicating whether the files are identical.
    """
    return filecmp.cmp(file1, file2)



class TestKratos(TestCase):

    root_folder_path = r'C:\checkouts\test_directory'
    subfolders_list = detect_subfolders(root_folder_path)

    def test_kratos_cases(self, subfolders_list=subfolders_list):
        """
        Run the test in every subfolder of the specified root folder.

        Parameters:
        - subfolders_list: A list of subfolders to run the test in.
        """

        for test_folder in subfolders_list:
            current_directory = os.getcwd()
            os.chdir(test_folder)
            name = test_folder.split('\\')[-1]
            with self.subTest(msg=name):
                run_test_in_folder()
                self.assertEqual(True, True)
            os.chdir(current_directory)
