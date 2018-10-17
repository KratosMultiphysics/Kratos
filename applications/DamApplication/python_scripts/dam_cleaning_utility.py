from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
import os

def CleanPreviousFiles(problem_path):

	file_ending_type = ".post.bin"
	CleanPreviousFileType(file_ending_type,problem_path)

	file_ending_type = ".post.res"
	CleanPreviousFileType(file_ending_type,problem_path)

	file_ending_type = ".post.msh"
	CleanPreviousFileType(file_ending_type,problem_path)

	file_ending_type = ".post.lst"
	CleanPreviousFileType(file_ending_type,problem_path)

	file_ending_type = "info.time"
	CleanPreviousFileType(file_ending_type,problem_path)


def CleanPreviousFileType(file_ending_type,problem_path):

    [os.remove(f) for f in os.listdir(problem_path) if f.endswith(file_ending_type)]
