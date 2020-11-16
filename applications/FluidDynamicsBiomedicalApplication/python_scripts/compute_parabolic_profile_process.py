# Importing the Kratos Library
import KratosMultiphysics

# Import applications
import KratosMultiphysics.FluidDynamicsBiomedicalApplication as FluidDynamicsBiomedicalApplication

# # Import base class file
# from FluidDynamicsBiomedicalApplication import WssStatisticsUtilities as wss_stadistics


def Factory(settings, model):
    if(type(settings) != KratosMultiphysics.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return ComputeParabolicProfileProcess(model, settings["Parameters"])

class ComputeParabolicProfileProcess(KratosMultiphysics.Process):
    """
    Auxiliary base class to compute parabolic inlet profile
    A derived class needs to be implemented to be able to use this functionality, as calling the base class alone is not enough.
    """

    def __init__(self, model, settings):
        # Base class constructor call
        KratosMultiphysics.Process.__init__(self)

        # Check default settings
        default_settings = KratosMultiphysics.Parameters("""
        {
            "inlet_model_part": "",
            "calculate_parabolic_profile": true,
        }
        """)
        settings.ValidateAndAssignDefaults(default_settings)

        # Save model and settings containers
        self.model = model
        self.settings = settings

    def ExecuteInitialize(self):
        # Compute the normal on the nodes of interest
        # Note that this overwrites the existent nodal normal values
        # Also note that if there is a slip condition that shares part of the WSS model part the corner normals could be altered
        # TODO: Improve the NormalCalculationUtils to accept alternative storage variables (to be discussed)
        inlet_model_part = self.model.GetModelPart(self.settings["inlet_model_part"].GetString())
        KratosMultiphysics.NormalCalculationUtils().CalculateOnSimplex(
            inlet_model_part,
            inlet_model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE])
        inlet_value=self.model.GetModelPart(self.settings["modulus"].GetString())  # TODO: inlet value is time dependent....
        if (self.settings["calculate_parabolic_profile"].GetBool()):
            FluidDynamicsBiomedicalApplication.CalculateParabolicProfile.ParabolicProfileMain(inlet_model_part,inlet_value)
