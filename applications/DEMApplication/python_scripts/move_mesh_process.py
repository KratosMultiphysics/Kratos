from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
# Importing the Kratos Library
import KratosMultiphysics
import KratosMultiphysics.DEMApplication as KratosDem
from KratosMultiphysics.DEMApplication import *


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

        self.model_part.SetValue(KratosMultiphysics.DEMApplication.LINEAR_VELOCITY, [0.0,0.0,0.0])
        self.model_part.SetValue(KratosMultiphysics.VELOCITY_PERIOD, 0.0)
        self.model_part.SetValue(KratosMultiphysics.ANGULAR_VELOCITY, [0.0,0.0,50.0])
        self.model_part.SetValue(KratosMultiphysics.ROTATION_CENTER, [0.0,0.0,0.0])
        self.model_part.SetValue(KratosMultiphysics.ANGULAR_VELOCITY_PERIOD, 0.0)
        self.model_part.SetValue(KratosMultiphysics.DEMApplication.VELOCITY_START_TIME, 0.0)
        self.model_part.SetValue(KratosMultiphysics.DEMApplication.VELOCITY_STOP_TIME, 1000.0)
        self.model_part.SetValue(KratosMultiphysics.DEMApplication.ANGULAR_VELOCITY_START_TIME, 0.0)
        self.model_part.SetValue(KratosMultiphysics.DEMApplication.ANGULAR_VELOCITY_STOP_TIME, 1000.0)
        self.model_part.SetValue(KratosMultiphysics.FIXED_MESH_OPTION, False)
        self.model_part.SetValue(KratosMultiphysics.DEMApplication.RIGID_BODY_MOTION, True)
        self.model_part.SetValue(KratosMultiphysics.DEMApplication.FREE_BODY_MOTION, False)
        self.model_part.SetValue(KratosMultiphysics.DEMApplication.IS_GHOST, False)
        self.model_part.SetValue(KratosMultiphysics.DEMApplication.TOP, False)
        self.model_part.SetValue(KratosMultiphysics.DEMApplication.BOTTOM, False)
        self.model_part.SetValue(KratosMultiphysics.DEMApplication.FORCE_INTEGRATION_GROUP, 0)

        self.mesh_motion = KratosMultiphysics.DEMApplication.DEMFEMUtilities()

    def ExecuteFinalizeSolutionStep(self):
        time = self.model_part.ProcessInfo[KratosMultiphysics.TIME]

        print("time: " + str(time))
        print("self.model_part: " + str(self.model_part.NumberOfNodes()))

        self.mesh_motion.MoveAllMeshes(self.model_part, time, self.time_step)

        print(str(self.model_part.GetNode(4).X) + " " + str(self.model_part.GetNode(4).Y) + " " + str(self.model_part.GetNode(4).Z))