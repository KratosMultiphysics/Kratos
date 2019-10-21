from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
# Importing the Kratos Library
import KratosMultiphysics
import KratosMultiphysics.DEMApplication as KratosDem


def Factory(settings, model):
    if(type(settings) != KratosMultiphysics.Parameters):
        raise Exception("Expected input shall be a Parameters object, encapsulating a json string")
    return MoveMeshProcess(model, settings["Parameters"])

class MoveMeshProcess(KratosMultiphysics.Process):

    def __init__(self, model, params):
        KratosMultiphysics.Process.__init__(self)

        ## Settings string in json format
        default_parameters = KratosMultiphysics.Parameters("""{
            "model_part_name"                : "",
            "time_step"                      : 0.1
        }""")

        ## Overwrite the default settings with user-provided parameters
        self.params = params
        self.params.ValidateAndAssignDefaults(default_parameters)

        self.model_part = model[self.params["model_part_name"].GetString()]

        self.time_step = self.params["time_step"].GetDouble()

        self.model_part.SetValue(RIGID_BODY_MOTION, True)
        self.model_part.SetValue(FREE_BODY_MOTION, False)
        self.model_part.SetValue(FIXED_MESH_OPTION, False)
        self.model_part.SetValue(LINEAR_VELOCITY, [0.0,0.0,0.0])
        self.model_part.SetValue(VELOCITY_PERIOD, 0.0)
        self.model_part.SetValue(ANGULAR_VELOCITY, [0.0,0.0,50.0])
        self.model_part.SetValue(ROTATION_CENTER, [0.0,0.0,0.0])
        self.model_part.SetValue(ANGULAR_VELOCITY_PERIOD, 0.0)
        self.model_part.SetValue(VELOCITY_START_TIME, 0.0)
        self.model_part.SetValue(VELOCITY_STOP_TIME, 1000.0)
        self.model_part.SetValue(ANGULAR_VELOCITY_START_TIME, 0.0)
        self.model_part.SetValue(ANGULAR_VELOCITY_STOP_TIME, 1000.0)
        self.model_part.SetValue(IS_GHOST, False)
        self.model_part.SetValue(TOP, False)
        self.model_part.SetValue(BOTTOM, False)
        self.model_part.SetValue(FORCE_INTEGRATION_GROUP, 0)

        self.mesh_motion = DEMFEMUtilities()

    def ExecuteFinalizeSolutionStep(self):
        time = self.model_part.ProcessInfo[KratosMultiphysics.TIME]

        center_position = [0,0,0]
        linear_velocity_changed = [0,0,0]
        angular_velocity_changed = [0,0,0]
        new_axes1 = [0,0,0]
        new_axes2 = [0,0,0]
        new_axes3 = [0,0,0]

        self.mesh_motion.MoveAllMeshes(self.model_part, time, self.time_step)


        print(self.model_part.GetNode(4))