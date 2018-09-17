from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7

# Importing the Kratos Library
import KratosMultiphysics


def Factory(settings, Model):
    if(type(settings) != KratosMultiphysics.Parameters):
        raise Exception("Expected input shall be a Parameters object, encapsulating a json string")
    return ALEBoundaryConditionProcess(Model, settings["Parameters"])

class ALEBoundaryConditionProcess(KratosMultiphysics.Process):
    def __init__(self, Model, settings):

        default_settings = KratosMultiphysics.Parameters(
            """
            {
                "help"            : "This process imposes the ALE boundary condition by fixing the MESH_VELOCITY",
                "model_part_name" : "please_specify_model_part_name",
                "constrained"     : [true,true,true],
            }
            """
        )

        settings.ValidateAndAssignDefaults(default_settings)
        settings.RemoveValue("help")

        KratosMultiphysics.Process.__init__(self)
        self.model_part = Model[settings["model_part_name"].GetString()]
        self.settings = settings

    def ExecuteInitialize(self):
        if not self.settings["constrained"].IsArray():
            raise Exception('"constrained" has to be be provided as an Array!')
        if self.settings["constrained"].size() != 3:
            raise Exception('"constrained" requires three entries!')

        mesh_vel_components = [
            KratosMultiphysics.MESH_VELOCITY_X,
            KratosMultiphysics.MESH_VELOCITY_Y,
            KratosMultiphysics.MESH_VELOCITY_Z
        ]

        for i in range(3):
            is_constrained = self.settings["constrained"][i].GetBool()
            KratosMultiphysics.VariableUtils().ApplyFixity(mesh_vel_components[i],
                                                           is_constrained,
                                                           self.model_part.Nodes)

    def ExecuteBeforeSolutionLoop(self):
        pass

    def ExecuteInitializeSolutionStep(self):
        pass

    def ExecuteFinalizeSolutionStep(self):
        pass

    def ExecuteBeforeOutputStep(self):
        pass

    def ExecuteAfterOutputStep(self):
        pass

    def ExecuteFinalize(self):
        pass
