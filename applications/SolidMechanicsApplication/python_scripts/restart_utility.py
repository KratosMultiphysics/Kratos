from __future__ import unicode_literals, print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
# importing the Kratos Library
from KratosMultiphysics import *
from KratosMultiphysics.SolidMechanicsApplication import *
CheckForPreviousImport()
import os


class RestartUtility:
    #

    def __init__(self, model_part, problem_path, problem_name):

        self.model_part = model_part

        # set restart flags
        self.load_restart_flag = False
        self.save_restart_flag = False

        # set problem name
        self.problem_name = problem_name

        # set problem path
        self.problem_path = problem_path

        # set serializer flag
        self.serializer_flag = "SERIALIZER_NO_TRACE"
        # self.serializer_flag = "SERIALIZER_TRACE_ERROR"

    #
    def Load(self, restart_step):

        self.load_restart_flag = True
        restart_path = os.path.join(self.problem_path, self.problem_name + "_" + str(restart_step))

        if(os.path.exists(restart_path) == False):
            print(" RESTART file do not exists , check the RestartStep selected ")

        kratos_serializer_variable = globals()[self.serializer_flag]
        serializer = Serializer(restart_path, kratos_serializer_variable)

        serializer.Load("ModelPart", self.model_part)

        print(" RESTART LOADED : [STEP: ", restart_step, "]")

    #
    def CleanPreviousFileType(self, file_ending_type):

        filelist = [f for f in os.listdir(self.problem_path) if f.endswith(file_ending_type)]

        for f in filelist:
            try:
                os.remove(f)
            except WindowsError:
                pass

    #
    def CleanPosteriorFileType(self, restart_step, file_ending_type):

        total_files = 0
        filelist = []
        for f in os.listdir(self.problem_path):
            if(f.endswith(file_ending_type)):
                total_files = total_files + 1

        total_files = total_files + 100  # arbitrary number to ensure to remove all files

        for rfile in range(0, total_files):
            for f in os.listdir(self.problem_path):
                step_time = restart_step + rfile
                if(f.endswith("_" + str(step_time) + file_ending_type)):
                    filelist.append(f)

        for f in filelist:
            try:
                os.remove(f)
            except WindowsError:
                pass

    #
    def CleanPreviousFiles(self):

        print("Start: -remove previous problem files-")

        # remove previous results:
        file_ending_type = ".post.bin"
        self.CleanPreviousFileType(file_ending_type)

        file_ending_type = ".post.res"
        self.CleanPreviousFileType(file_ending_type)

        file_ending_type = ".post.msh"
        self.CleanPreviousFileType(file_ending_type)

        # remove previous graph files:
        file_ending_type = ".png"
        self.CleanPreviousFileType(file_ending_type)

        # remove previous restart files:
        file_ending_type = ".rest"
        self.CleanPreviousFileType(file_ending_type)

    #
    def CleanPosteriorFiles(self, restart_step):

        print("Restart: -remove restart step posterior problem files-")

        # set load restart step
        self.restart_step = restart_step + 1

        # remove previous results after restart:
        file_ending_type = ".post.bin"
        self.CleanPosteriorFileType(self.restart_step, file_ending_type)

        file_ending_type = ".post.res"
        self.CleanPosteriorFileType(self.restart_step, file_ending_type)

        file_ending_type = ".post.msh"
        self.CleanPosteriorFileType(self.restart_step, file_ending_type)

        # remove previous graphs after restart:

        file_ending_type = ".png"
        self.CleanPosteriorFileType(self.restart_step, file_ending_type)

        # remove previous restart files:
        file_ending_type = ".rest"
        self.CleanPosteriorFileType(self.restart_step, file_ending_type)

    #
    def Save(self, current_time, current_step, current_id):

        label = current_step + 1

        restart_path = os.path.join(self.problem_path, self.problem_name + "_" + str(label))

        kratos_serializer_variable = globals()[self.serializer_flag]
        serializer = Serializer(restart_path, kratos_serializer_variable)

        serializer.Save("ModelPart", self.model_part)
        print(" RESTART PRINTED: [ID: ", label, "] [STEP: ", current_step, "] [TIME: ", current_time, "]")

    #
