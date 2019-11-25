# Importing the Kratos Library
import KratosMultiphysics

# Import applications
import KratosMultiphysics.FluidDynamicsBiomedicalApplication as FluidDynamicsBiomedicalApplication

# # Import base class file
# from FluidDynamicsBiomedicalApplication import WssStatisticsUtilities as wss_stadistics


def Factory(settings, model):
    if(type(settings) != KratosMultiphysics.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return ComputeWssStatisticsProcess(model, settings["Parameters"])


class ComputeWssStatisticsProcess(KratosMultiphysics.Process):
    """
    Auxiliary base class to output total WSS forces
    A derived class needs to be implemented to be able to use this functionality, as calling the base class alone is not enough.
    """

    def __init__(self, model, settings):
        # Base class constructor call
        KratosMultiphysics.Process.__init__(self)

        # Check default settings
        default_settings = KratosMultiphysics.Parameters("""
        {
            "model_part_name": "please_specify_skin_model_part_name",
            "calculate_wss": true,
            "calculate_osi": true
        }
        """)
        settings.ValidateAndAssignDefaults(default_settings)

        # Save model and settings containers
        self.model = model
        self.settings = settings

    def ExecuteFinalizeSolutionStep(self):
        model_part = self.model.GetModelPart(self.settings["model_part_name"].GetString())

        if (self.settings["calculate_wss"].GetBool()):
            FluidDynamicsBiomedicalApplication.WssStatisticsUtilities.CalculateWSS(model_part)

    def ExecuteFinalize(self):
        model_part = self.model.GetModelPart(self.settings["model_part_name"].GetString())

        if (self.settings["calculate_osi"].GetBool()):
            FluidDynamicsBiomedicalApplication.WssStatisticsUtilities.CalculateOSI(model_part)
            FluidDynamicsBiomedicalApplication.WssStatisticsUtilities.CalculateTWSS(model_part)