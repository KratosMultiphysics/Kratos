from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7

import os
#
import ListParameters  as list_variables

# set problem name
problem_name = list_variables.problem_name

# set problem path
problem_path = list_variables.problem_path

# set print list
if(list_variables.print_lists == "True"):
    print_lists = True
else:
    print_lists = False

# set output mode
if(list_variables.output_mode == "Binary"):
    output_mode = ".post.bin"
if(list_variables.output_mode == "Ascii"):
    output_mode = ".post.res"

# initialize lists
file_list = list_variables.file_list

header_in_list = []
listprint = []


# ----------------------------------------------------------------#
#
def RemoveListFiles():

    # remove previous list files:
    if(os.path.exists(problem_path) == False):
        print(" Problem Path do not exists , check the Problem Path selected ")
    else:
        filelist = [f for f in os.listdir(problem_path) if f.endswith(".lst")]
        for f in filelist:
            if(os.path.exists(f)):
                try:
                    os.remove(f)
                except OSError:
                    pass

# ----------------------------------------------------------------#
#
def FileExists(file_name):

    file_exists = False

    if(os.path.exists(file_name)):
        file_exists = True

    return file_exists

# ----------------------------------------------------------------#
#
def ReBuildListFiles():

    RemoveListFiles()

    # Rebuild List Files from existing problem build
    if(print_lists):

        file_id   = []

        for f in os.listdir(problem_path):

            if(f.endswith(output_mode)):

                #if f.name = problem_tested_145.post.bin
                file_parts = f.split('_')  # you get ["problem","tested","145.post.bin"]
                num_parts  = len(file_parts)

                end_parts  = file_parts[num_parts-1].split(".") # you get ["145","post","bin"]
                print_id   = end_parts[0] # you get "145"

                file_id.append(int(print_id))


        file_id.sort()

        num_list_files = len(file_id)

        for lfile in range(0, num_list_files):

            print_id   = file_id[lfile]

            num_list_files = len(file_list)

            for lfile in range(0, num_list_files):
                if(file_list[lfile] == listprint[lfile]):

                    problempath = os.path.join(problem_path, problem_name + "_" + str(file_list[lfile]) + ".post.lst")

                    if(FileExists(problempath) == False):
                        header_in_list[lfile]

                    listfile = open(problempath, "a")

                    if(header_in_list[lfile]):
                        problemname = "Multiple\n"
                        listfile.write(problemname)
                        header_in_list[lfile] = False

                    problemname = problem_name + "_" + str(print_id) + output_mode + "\n"
                    listfile.write(problemname)
                    listfile.close()
                    listprint[lfile] = 1
                else:
                    listprint[lfile] = listprint[lfile] + 1

# ----------------------------------------------------------------#




# check if a list with all files is required
all_files_list = False
num_list_files = len(file_list)

for lfile in range(0, num_list_files):
    if(file_list[lfile] == 1):
        all_files_list = True

# add a list with all files:
if(all_files_list == False):
    file_list.append(1)

for lfile in range(0, num_list_files):
    file_list.append(file_list[lfile])

for lfile in file_list:
    listprint.append(1)
    header_in_list.append(True)


ReBuildListFiles()
