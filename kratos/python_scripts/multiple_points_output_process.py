from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7

# Importing the Kratos Library
import KratosMultiphysics

# other imports
from point_output_process import PointOutputProcess

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
        KratosMultiphysics.Process.__init__(self)

        default_settings = KratosMultiphysics.Parameters('''{
            "help"              : "This process writes several points to a file. Internally it holds objects of type PointOutputProcess",
            "model_part_name"   : "",
            "entity_type"       : "element",
            "positions"         : [[]],
            "output_variables"  : [],
            "historical_value"  : true,
            "print_format"      : "",
            "output_file_settings": {}
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

        if params["output_file_settings"]["file_name"].GetString().endswith(".dat"):
            params["output_file_settings"]["file_name"].SetString(params["output_file_settings"]["file_name"].GetString()[:-4])

        # Create the individual point_output_processes
        for i in range(num_points):
            point_proc_params = params.Clone()

            for j in range(3):
                position_vec[j] = positions[i,j]
            point_proc_params["position"].SetVector(position_vec)
            point_proc_params["output_file_settings"]["file_name"].SetString(params["output_file_settings"]["file_name"].GetString() + "_" + str(i+1))

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
