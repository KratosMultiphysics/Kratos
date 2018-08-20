from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7

# Importing the Kratos Library
import KratosMultiphysics

# other imports
from multiple_points_output_process import MultiplePointsOutputProcess

def Factory(settings, Model):
    if(type(settings) != KratosMultiphysics.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return LineOutputProcess(Model, settings["Parameters"])

class LineOutputProcess(KratosMultiphysics.Process):
    """This process writes output for several points along a line to a file
    Internally it holds an object of type "MultiplePointsOutputProcess"
    Usage:
        - specifying start and end point defining the line
        - specifying the number of sampling points along the line (start and end points will be included)
    """
    def __init__(self, model, params):

        default_settings = KratosMultiphysics.Parameters('''{
            "help"              : "This process writes output for several points along a line to a file. Internally it holds an object of type MultiplePointsOutputProcess",
            "model_part_name"   : "",
            "entity_type"       : "element",
            "start_point"       : [],
            "end_point"         : [],
            "sampling_points"   : 3,
            "output_variables"  : [],
            "print_format"      : "",
            "output_file_settings": {}
        }''')

        params.ValidateAndAssignDefaults(default_settings)

        # retrieving the positions defining the line entity as well as the number of sampling points
        start_point_position = params["start_point"].GetVector()
        if start_point_position.Size() != 3:
            raise Exception('The start point position has to be provided with 3 coordinates!')
        end_point_position = params["end_point"].GetVector()
        if end_point_position.Size() != 3:
            raise Exception('The end point position has to be provided with 3 coordinates!')
        number_of_sampling_points = params["sampling_points"].GetInt()

        # check the entity type
        entity_type = params["entity_type"].GetString()
        if not entity_type == "element" or entity_type == "condition":
            raise Exception('"entity_type" can be either "element" or "condition"!')

        if number_of_sampling_points <= 0:
            raise Exception('The number of sampling points has to be larger than 0!')

        elif number_of_sampling_points == 1:
            positions = KratosMultiphysics.Matrix(1, 3)
            positions[0,0] = (start_point_position[0] + end_point_position[0] ) / 2.
            positions[0,1] = (start_point_position[1] + end_point_position[1] ) / 2.
            positions[0,2] = (start_point_position[2] + end_point_position[2] ) / 2.

        else:
            # setup the parametric space for the internal points on the line
            lower_bound = 0.0
            upper_bound = 1.0
            parametrized_internal_points = [lower_bound + x*(upper_bound-lower_bound)/(number_of_sampling_points-1) for x in range(number_of_sampling_points)]

            # determining the positions of the output points
            direction_vector = [x - y for x, y in zip(end_point_position, start_point_position)]

            positions = KratosMultiphysics.Matrix(len(parametrized_internal_points), 3)
            for k in range(len(parametrized_internal_points)):
                current_position = [x + parametrized_internal_points[k]*y for x, y in zip(start_point_position, direction_vector)]

                positions[k,0] = current_position[0]
                positions[k,1] = current_position[1]
                positions[k,2] = current_position[2]

        params.RemoveValue("start_point")
        params.RemoveValue("end_point")
        params.RemoveValue("sampling_points")

        params.AddEmptyValue("positions")
        params["positions"].SetMatrix(positions)

        self.multiple_points_output_process = MultiplePointsOutputProcess(model, params)

    def ExecuteInitialize(self):
        self.multiple_points_output_process.ExecuteInitialize()

    def ExecuteBeforeSolutionLoop(self):
        self.multiple_points_output_process.ExecuteBeforeSolutionLoop()

    def ExecuteInitializeSolutionStep(self):
        self.multiple_points_output_process.ExecuteInitializeSolutionStep()

    def ExecuteFinalizeSolutionStep(self):
        self.multiple_points_output_process.ExecuteFinalizeSolutionStep()

    def ExecuteBeforeOutputStep(self):
        self.multiple_points_output_process.ExecuteBeforeOutputStep()

    def ExecuteAfterOutputStep(self):
        self.multiple_points_output_process.ExecuteAfterOutputStep()

    def ExecuteFinalize(self):
        self.multiple_points_output_process.ExecuteFinalize()
