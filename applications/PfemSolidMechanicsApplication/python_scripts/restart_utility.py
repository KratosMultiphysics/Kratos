from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
# importing the Kratos Library
from KratosMultiphysics import *
from KratosMultiphysics.SolidMechanicsApplication import *
CheckForPreviousImport()
import os


class RestartUtility:
    #

    def __init__(self, model_part, problem_path, problem_name):


        self.model_part =  model_part

        # set restart flags
        self.load_restart_flag = False
        self.save_restart_flag = False

        # set problem name
        self.problem_name = problem_name

        # set problem path
        self.problem_path = problem_path

        print(" problem path ", self.problem_path)

        # set serializer flag
        self.serializer_flag = "SERIALIZER_NO_TRACE"      # binary
        # self.serializer_flag = "SERIALIZER_TRACE_ERROR" # ascii
        # self.serializer_flag = "SERIALIZER_TRACE_ALL"   # ascii

    #
    def Load(self, restart_step):

        self.load_restart_flag = True
        restart_path = os.path.join(self.problem_path, self.problem_name + "_" + str(restart_step))

        if(os.path.exists(restart_path+".rest") == False):
            print("::Restart Utility:: RESTART file does not exist , check the RestartStep selected ")

        kratos_serializer_variable = globals()[self.serializer_flag]
        serializer = Serializer(restart_path, kratos_serializer_variable)

        serializer.Load("ModelPart", self.model_part)

        print(" ::[Restart Utility]:: RESTART STEP",restart_step,"LOADED]" )

    #
    def CleanPreviousFileType(self, file_ending_type):

        if(os.path.exists(self.problem_path) == False):
            print("::Restart Utility:: Problem Path does not exist , check the Problem Path selected ")
        else:
            filelist = [f for f in os.listdir(self.problem_path) if f.endswith(file_ending_type)]

            for f in filelist:
                try:
                    os.remove(f)
                except WindowsError:
                    pass

    #
    def CleanPosteriorFileType(self, restart_step, file_ending_type):

        if(os.path.exists(self.problem_path) == False):
            print(" ::[Restart Utility]:: Problem Path does not exist , check the Problem Path selected ")
        else:
            filelist = []
            for f in os.listdir(self.problem_path):
                if(f.endswith(file_ending_type)):
                    #if f.name = problem_tested_145.post.bin
                    file_parts = f.split('_')  # you get ["problem","tested","145.post.bin"]
                    num_parts  = len(file_parts) 
                    
                    end_parts  = file_parts[num_parts-1].split(".") # you get ["145","post","bin"]
                    print_id   = end_parts[0] # you get "145"

                    if( int(print_id)>restart_step ):
                        filelist.append(f)

            for f in filelist:
                try:
                    os.remove(f)
                except WindowsError:
                    pass

    #
    def CleanPreviousFiles(self):

        #print(" ::[Restart Utility]:: Remove Previous Files")

        # remove previous results:
        file_ending_type = ".post.bin"
        self.CleanPreviousFileType(file_ending_type)

        file_ending_type = ".post.res"
        self.CleanPreviousFileType(file_ending_type)

        file_ending_type = ".post.msh"
        self.CleanPreviousFileType(file_ending_type)

        # remove previous graph files:
        file_ending_type = ".graph.png"
        self.CleanPreviousFileType(file_ending_type)

        # remove previous data files:
        file_ending_type = ".post.csv"
        self.CleanPreviousFileType(file_ending_type)

        # remove previous restart files:
        file_ending_type = ".rest"
        self.CleanPreviousFileType(file_ending_type)

    #
    def CleanPosteriorFiles(self, restart_step):

        #print(" ::[Restart Utility]:: Clean Post Restart Files)")

        # remove posterior results after restart:
        file_ending_type = ".post.bin"
        self.CleanPosteriorFileType(restart_step, file_ending_type)

        file_ending_type = ".post.res"
        self.CleanPosteriorFileType(restart_step, file_ending_type)

        file_ending_type = ".post.msh"
        self.CleanPosteriorFileType(restart_step, file_ending_type)

        # remove posterior graphs after restart:

        file_ending_type = "graph.png"
        self.CleanPosteriorFileType(restart_step+1, file_ending_type)

        # remove posterior restart files:
        file_ending_type = ".rest"
        self.CleanPosteriorFileType(restart_step+1, file_ending_type)

    #
    def Save(self, current_time, current_step, current_id):

        label = current_id

        restart_path = os.path.join(self.problem_path, self.problem_name + "_" + str(label))

        kratos_serializer_variable = globals()[self.serializer_flag]
        serializer = Serializer(restart_path, kratos_serializer_variable)

        serializer.Save("ModelPart", self.model_part)
        #print(" ::[Restart Utility]:: RESTART PRINTED: [ID: ", label, "] [STEP: ", current_step, "] [TIME: ", current_time, "]")

    #
