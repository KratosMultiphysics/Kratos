from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7

# Importing the Kratos Library
import KratosMultiphysics

# other imports
from point_output_process import PointOutputProcess
import os

def Factory(settings, Model):
    if(type(settings) != KratosMultiphysics.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return MultiplePointsOutputProcess(Model, settings["Parameters"])

class MultiplePointsOutputProcess(KratosMultiphysics.Process):
    """This process writes several points to a file
    Internally it holds objects of type "PointOutputProcess"
    Usage:
        - directly by specifying the locations of the points for which output is wanted
        - inside other processes that create points for which output is required
          (e.g. LinePointsOutputProcess)
    """
    def __init__(self, model, params):

        default_settings = KratosMultiphysics.Parameters('''{
            "positions"          : [[]],
            "model_part_name"    : "",
            "output_file_name"   : "",
            "output_variables"   : [],
            "entity_type"        : "element",
            "save_output_file_in_folder"  : true,
            "output_folder_relative_path" : "TabularResults"
        }''')

        params.ValidateAndAssignDefaults(default_settings)

        positions = params["positions"].GetMatrix()

        num_points = positions.Size1()
        if num_points == 0:
            raise Exception('No positions were specified"')
        if positions.Size2() != 3:
            raise Exception('The positions have to be provided with 3 coordinates!')

        self.point_output_processes = []
        params.RemoveValue("positions")
        params.AddEmptyValue("position")
        position_vec = KratosMultiphysics.Vector(3)

        # setting up the output_file
        raw_path, output_file_name_base = os.path.split(params["output_file_name"].GetString())
        if output_file_name_base == "":
            raise Exception('No "output_file_name" was specified!')

        if params["save_output_file_in_folder"].GetBool():
            if params["output_folder_relative_path"].GetString() == "":
                raise Exception('No "save_output_file_in_folder" was specified!')
            else:
                output_folder_relative_path = os.path.join(params["output_folder_relative_path"].GetString(),
                                                            "mp_" + output_file_name_base)
                if raw_path != "":
                    warn_msg  = 'Relative path "'+ raw_path +'" contained wrongly in "output_file_name" : "'+ params["output_file_name"].GetString() +'"\n'
                    warn_msg += 'Use "output_folder_relative_path" to specify correctly\n'
                    warn_msg += 'Using the default relative path "' + output_folder_relative_path + '" instead'
                    KratosMultiphysics.Logger.PrintWarning("MultiplePointsOutputProcess", warn_msg)
        else:
            if raw_path != "":
                warn_msg  = 'Relative path "'+ raw_path +'" contained wrongly in "output_file_name": "'+ params["output_file_name"].GetString() +'"\n'
                warn_msg += 'Use the "save_output_file_in_folder" and "output_folder_relative_path" to specify correctly\n'
                warn_msg += 'Using the current directory instead'
                KratosMultiphysics.Logger.PrintWarning("MultiplePointsOutputProcess", warn_msg)

        # Create the individual point_output_processes
        for i in range(num_points):
            point_proc_params = params.Clone()

            for j in range(3):
                position_vec[j] = positions[i,j]
            point_proc_params["position"].SetVector(position_vec)

            output_file_name = output_file_name_base + "_" + str(i+1)

            point_proc_params["output_file_name"].SetString(output_file_name)
            if params["save_output_file_in_folder"].GetBool():
                point_proc_params["output_folder_relative_path"].SetString(output_folder_relative_path)
            self.point_output_processes.append(PointOutputProcess(model, point_proc_params))

    def ExecuteInitialize(self):
        for proc in self.point_output_processes:
            proc.ExecuteInitialize()

    def ExecuteBeforeSolutionLoop(self):
        for proc in self.point_output_processes:
            proc.ExecuteBeforeSolutionLoop()

    def ExecuteInitializeSolutionStep(self):
        for proc in self.point_output_processes:
            proc.ExecuteInitializeSolutionStep()

    def ExecuteFinalizeSolutionStep(self):
        for proc in self.point_output_processes:
            proc.ExecuteFinalizeSolutionStep()

    def ExecuteBeforeOutputStep(self):
        for proc in self.point_output_processes:
            proc.ExecuteBeforeOutputStep()

    def ExecuteAfterOutputStep(self):
        for proc in self.point_output_processes:
            proc.ExecuteAfterOutputStep()

    def ExecuteFinalize(self):
        for proc in self.point_output_processes:
            proc.ExecuteFinalize()
