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

    for f in os.listdir(problem_path):
        if f.endswith(file_ending_type):
            os.remove(f)
