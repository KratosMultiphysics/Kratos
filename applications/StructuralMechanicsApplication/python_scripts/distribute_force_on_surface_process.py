import KratosMultiphysics as KM
import KratosMultiphysics.StructuralMechanicsApplication as KSM


def Factory(settings, Model):
    if not isinstance(settings, KM.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return DistributeForceOnSurfaceProcess(Model, settings["Parameters"])


class DistributeForceOnSurfaceProcess(KM.Process):
    """This process distributes a force on surface load conditions belonging to a submodelpart.
    The force is distributed according to the surface area.
    """
    def __init__(self, model, settings):
        """ The default constructor of the class

        Arguments:
        model -- the container of the different model parts.
        settings -- Kratos parameters containing process settings.
        """
        KM.Process.__init__(self)

        # The value can be a double or a string (function)
        default_settings = KM.Parameters("""
        {
            "help"                 : "This process distributes a force on surface load conditions belonging to a submodelpart. The force is distributed according to the surface area.",
            "model_part_name"      : "please_specify_model_part_name",
            "interval"             : [0.0, 1e30],
            "modulus"              : 1.0,
            "direction"            : [1.0, 0.0, 0.0]
        }
        """)

        # Assign this here since it will change the "interval" prior to validation
        self.interval = KM.IntervalUtility(settings)

        settings.ValidateAndAssignDefaults(default_settings)

        self.variable_utils = KM.VariableUtils()

        self.model_part = model.GetModelPart(settings["model_part_name"].GetString())
        modulus = settings["modulus"].GetDouble()
        self.force = [
            settings["direction"][0].GetDouble() * modulus,
            settings["direction"][1].GetDouble() * modulus,
            settings["direction"][2].GetDouble() * modulus
        ]

    def ExecuteInitializeSolutionStep(self):
        """Calculate the total surface area of the conditions and assign the SURFACE_LOAD to the conditions
        """
        current_time = self.model_part.ProcessInfo[KM.TIME]

        if self.interval.IsInInterval(current_time):

            total_area = sum(x.GetGeometry().Area() for x in self.model_part.Conditions)

            force_by_area = KM.Vector(3)
            force_by_area[0] = self.force[0] / total_area
            force_by_area[1] = self.force[1] / total_area
            force_by_area[2] = self.force[2] / total_area

            for condition in self.model_part.Conditions:
                area = condition.GetGeometry().Area()
                condition.SetValue(KSM.SURFACE_LOAD, force_by_area * area)
