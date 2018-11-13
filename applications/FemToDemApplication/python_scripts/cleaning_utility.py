from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
import os


class CleaningUtility:

    def __init__(self, problem_path):

        self.problem_path = problem_path


    def CleanPreviousFiles(self):

        print("Removing previous problem files")

        file_ending_type = ".post.bin"
        self.CleanPreviousFileType(file_ending_type)

        file_ending_type = ".post.res"
        self.CleanPreviousFileType(file_ending_type)

        file_ending_type = ".post.msh"
        self.CleanPreviousFileType(file_ending_type)


    def CleanPreviousFileType(self, file_ending_type):

        if(os.path.exists(self.problem_path) == False):
            print("Problem Path do not exists, check the Problem Path selected")
        else:
            filelist = [f for f in os.listdir(self.problem_path) if f.endswith(file_ending_type)]

            for f in filelist:
                try:
                    os.remove(f)
                except WindowsError:
                    pass
