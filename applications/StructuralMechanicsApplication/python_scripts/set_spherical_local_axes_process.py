# Importing the Kratos Library
import KratosMultiphysics as KM
import KratosMultiphysics.StructuralMechanicsApplication as SMA
from KratosMultiphysics import Logger

def Factory(settings, Model):
    if not isinstance(settings, KM.Parameters):
        raise Exception("Expected input shall be a Parameters object, encapsulating a json string")
    return SetSphericalLocalAxesProcess(Model, settings["Parameters"])

class SetSphericalLocalAxesProcess(KM.Process):
    def __init__(self, Model, settings):
        KM.Process.__init__(self)
        self.settings = settings;
        self.model_part = Model[self.settings["model_part_name"].GetString()]

    def ExecuteInitialize(self):
        default_settings = KM.Parameters(
            """{
                "model_part_name"            : "set_model_part_name",
                "spherical_reference_axis"   : [0.0,0.0,1.0],
                "spherical_central_point"    : [0.0,0.0,0.0]
            }""");

        self.settings.ValidateAndAssignDefaults(default_settings)
        KM.Process.__init__(self)
        # Let's compute the local axes
        Logger.PrintInfo("SetSphericalLocalAxesProcess:: ","Setting the oriented local axes...")
        self.settings.RemoveValue("model_part_name")
        SMA.SetSphericalLocalAxesProcess(self.model_part, self.settings).ExecuteInitialize()




