from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7

# Importing the Kratos Library
import KratosMultiphysics

def Factory(settings, Model):
    if(type(settings) != KratosMultiphysics.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return MoveMeshProcess(Model, settings["Parameters"])

class MoveMeshProcess(KratosMultiphysics.Process):
    """This process writes output for several points along a line to a file
    Internally it holds an object of type "MultiplePointsOutputProcess"
    Usage:
        - specifying start and end point defining the line
        - specifying the number of sampling points along the line (start and end points will be included)
    """
    def __init__(self, model, params):
        KratosMultiphysics.Process.__init__(self)

        default_parameters = KratosMultiphysics.Parameters('''{
            "model_part_name"   : "",
            "3d_curve" : {
              "degree": 0,
              "knot_vector": [ ],
              "active_range": [ ],
              "control_points": [[ 1, [0,0]]]
              },
            "reference_parameter"       : 1.0,
            "current_parameter"         : 1.0
        }''')

        ## Overwrite the default settings with user-provided parameters
        self.params = params
        self.params.ValidateAndAssignDefaults(default_parameters)

        self.model_part = model[self.params["model_part_name"].GetString()]

        print("here")

        self.reference_location = [0, 0, 0]
        self.reference_direction = [0, 0, 0]

    #def ExecuteInitialize(self):
    #    self.multiple_points_output_process.ExecuteInitialize()

    #def ExecuteBeforeSolutionLoop(self):
    #    self.multiple_points_output_process.ExecuteBeforeSolutionLoop()

    def ExecuteInitializeSolutionStep(self):
        print("here2")

        #Matrix trans = ZeroMatrix(2,2)
        #matrix rotation = ZeroMatrix(2,2)

        #for node in self.model_part.Nodes():
        #    new_coords = node.initial_coords * irgen
        #    node.X() = new_coords[0]


    #def ExecuteFinalizeSolutionStep(self):
    #    self.multiple_points_output_process.ExecuteFinalizeSolutionStep()

    #def ExecuteBeforeOutputStep(self):
    #    self.multiple_points_output_process.ExecuteBeforeOutputStep()

    #def ExecuteAfterOutputStep(self):
    #    self.multiple_points_output_process.ExecuteAfterOutputStep()

    #def ExecuteFinalize(self):
    #    self.multiple_points_output_process.ExecuteFinalize()
