# Importing the Kratos Library
import KratosMultiphysics as KM
import KratosMultiphysics.StructuralMechanicsApplication as SMA

def Factory(settings, Model):
    if not isinstance(settings, KM.Parameters):
        raise Exception("Expected input shall be a Parameters object, encapsulating a json string")
    return SetLocalAxesProcess(Model, settings["Parameters"])

class SetLocalAxesProcess(KM.Process):
    def __init__(self, Model, settings):
        KM.Process.__init__(self)
        self.settings = settings;
        self.model_part = Model[self.settings["model_part_name"].GetString()]

    def ExecuteInitialize(self):
        default_settings = KM.Parameters(
            """
            {
                "model_part_name"               : "set_model_part_name",
                "local_axes_coordinate_system"  : "cartesian",
                "cartesian_local_axis"          : [[1.0,0.0,0.0],[0.0,1.0,0.0],[0.0,0.0,1.0]],
                "cylindrical_generatrix_axis"   : [0.0,0.0,1.0],
                "cylindrical_generatrix_point"  : [0.0,0.0,0.0],
                "spherical_reference_axis"      : [0.0,0.0,1.0],
                "spherical_central_point"       : [0.0,0.0,0.0]
            }""");

        self.settings.ValidateAndAssignDefaults(default_settings)

        KM.Process.__init__(self)

        # Let's compute the local axes
        self.settings.RemoveValue("model_part_name")
        if (self.settings["local_axes_coordinate_system"].GetString() == "cartesian"):
            self.settings.RemoveValue("cylindrical_generatrix_axis")
            self.settings.RemoveValue("cylindrical_generatrix_point")
            self.settings.RemoveValue("spherical_reference_axis")
            self.settings.RemoveValue("spherical_central_point")
            SMA.SetLocalAxesUtility().SetLocalAxisCartesianSystem(self.model_part, self.settings)
        elif (self.settings["local_axes_coordinate_system"].GetString() == "cylindrical"):
            self.settings.RemoveValue("cartesian_local_axis")
            self.settings.RemoveValue("spherical_reference_axis")
            self.settings.RemoveValue("spherical_central_point")
            SMA.SetLocalAxesUtility().SetLocalAxisCylindricalSystem(self.model_part, self.settings)
        elif (self.settings["local_axes_coordinate_system"].GetString() == "spherical"):
            self.settings.RemoveValue("cartesian_local_axis")
            self.settings.RemoveValue("cylindrical_generatrix_axis")
            self.settings.RemoveValue("cylindrical_generatrix_point")
            SMA.SetLocalAxesUtility().SetLocalAxisSphericalSystem(self.model_part, self.settings)
        else:
            raise Exception("local_axes_coordinate_system not supported...")
    # def ExecuteInitializeSolutionStep(self):
    #     pass





