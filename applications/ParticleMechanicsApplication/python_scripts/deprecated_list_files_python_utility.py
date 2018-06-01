from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
# Class to manage list files and result files
from KratosMultiphysics import *
CheckForPreviousImport()


class ListFilesUtility:
    #

    def __init__(self, problem_path, problem_name, print_lists, output_mode):

        # set problem name
        self.problem_name = problem_name

        # set problem path
        self.problem_path = problem_path

        # set print list
        if(print_lists == "True"):
            self.print_lists = True
        else:
            self.print_lists = False

        # set output mode
        if(output_mode == "Binary"):
            self.output_mode = ".post.bin"
        if(output_mode == "Ascii"):
            self.output_mode = ".post.res"

        # initialize lists
        self.file_list = []
        self.header_in_list = []
        self.listprint = []

    #
    def Initialize(self, file_list):

        # check if a list with all files is required
        all_files_list = False

        num_list_files = len(file_list)

        for lfile in range(0, num_list_files):
            if(file_list[lfile] == 1):
                all_files_list = True

        # add a list with all files:
        if(all_files_list == False):
            self.file_list.append(1)

        for lfile in range(0, num_list_files):
            self.file_list.append(file_list[lfile])

        for lfile in self.file_list:
            self.listprint.append(1)
            self.header_in_list.append(True)

    #
    def FileExists(self, file_name):

        file_exists = False

        if(os.path.exists(file_name)):
            file_exists = True

        return file_exists

    #
    def ReBuildListFiles(self):

        self.RemoveListFiles()
                
        # Rebuild List Files from existing problem build
        if(self.print_lists):

            file_id   = []

            for f in os.listdir(self.problem_path):

                if(f.endswith(self.output_mode)):
                    
                    #if f.name = problem_tested_145.post.bin
                    file_parts = f.split('_')  # you get ["problem","tested","145.post.bin"]
                    num_parts  = len(file_parts) 
                    
                    end_parts  = file_parts[num_parts-1].split(".") # you get ["145","post","bin"]
                    print_id   = end_parts[0] # you get "145"

                    file_id.append(int(print_id))

  
            file_id.sort()

            num_list_files = len(file_id)

            print("::[List Utilities]:: Rebuild Lists (Previous Files:", num_list_files,")")

            for lfile in range(0, num_list_files):
               
                print_id   = file_id[lfile]
                   
                num_list_files = len(self.file_list)

                for lfile in range(0, num_list_files):
                    if(self.file_list[lfile] == self.listprint[lfile]):

                        problempath = os.path.join(self.problem_path, "_List_" + str(self.file_list[lfile]) + "_" + self.problem_name + ".post.lst")

                        if(self.FileExists(problempath) == False):
                            self.header_in_list[lfile] = True
                                
                        listfile = open(problempath, "a")

                        if(self.header_in_list[lfile]):
                            problemname = "Multiple\n"
                            listfile.write(problemname)
                            self.header_in_list[lfile] = False

                        problemname = self.problem_name + "_" + str(print_id) + self.output_mode + "\n"
                        listfile.write(problemname)

                        listfile.close()

                        self.listprint[lfile] = 1
                    else:
                        self.listprint[lfile] = self.listprint[lfile] + 1



    #
    def RemoveListFiles(self):

        # remove previous list files:
        if(os.path.exists(self.problem_path) == False):
            print(" Problem Path do not exists , check the Problem Path selected ")
        else:
            filelist = [f for f in os.listdir(self.problem_path) if f.endswith(".lst")]
            for f in filelist:
                if(os.path.exists(f)):
                    try:
                        os.remove(f)
                    except WindowsError:
                        pass

    #
    def PrintListFiles(self, current_id):

        if(os.path.exists(self.problem_path) == False):
            print(" Problem Path do not exists , check the Problem Path selected ")
        else:
            # print list files:
            if(self.print_lists):

                num_list_files = len(self.file_list)

                for lfile in range(0, num_list_files):

                    if(self.file_list[lfile] == self.listprint[lfile]):
                        problempath = os.path.join(self.problem_path, "_List_" + str(self.file_list[lfile]) + "_" + self.problem_name + ".post.lst")
                        listfile = open(problempath, "a")

                        if(self.header_in_list[lfile]):
                            problemname = "Multiple\n"
                            listfile.write(problemname)
                            self.header_in_list[lfile] = False
                            problemname = self.problem_name + "_" + str(0) + self.output_mode + "\n"
                            listfile.write(problemname)
                        
                        if( current_id !=0 ):
                            problemname = self.problem_name + "_" + str(current_id) + self.output_mode + "\n"
                            listfile.write(problemname)

                        listfile.close()
                        self.listprint[lfile] = 1
                    else:
                        self.listprint[lfile] = self.listprint[lfile] + 1

    #

