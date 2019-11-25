# Importing the Kratos Library
import KratosMultiphysics

# Import applications
import KratosMultiphysics.FluidDynamicsBiomedicalApplication as FluidDynamicsBiomedicalApplication

# Import base class file
from FluidDynamicsBiomedicalApplication import WssStatisticsUtilities as wss_stadistics


def Factory(settings, model):
    if(type(settings) != KratosMultiphysics.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return ComputeWSSProcess(model, settings["Parameters"])


class ComputeWSSProcess(KratosMultiphysics.Process):
    """
    Auxiliary base class to output total WSS forces
    A derived class needs to be implemented to be able to use this functionality, as calling the base class alone is not enough.
    """
    def __init__(self, model, params ):
        """
        Auxiliary class to output total WSS forces.
        """
        KratosMultiphysics.Process.__init__(self)

    def __init__(self, model, settings):
        KratosMultiphysics.Process.__init__(self)
        default_settings = KratosMultiphysics.Parameters("""
            {
                "model_part_name": "please_specify_skin_model_part_name",
                "calculate_wss": true,
                "calculate_osi": true
            }
            """)

        model_part_name = self.params["model_part_name"].GetString()

    def ExecuteFinalizeSolutionStep(self):
        self.calculate_wss = self.params["calculate_wss"].GetBool()
        if (self.calculate_wss):
            wss_stadistics.CalculateWSS(model_part_name)

    def ExecuteFinalize(self):
        self.calculate_osi = self.params["calculate_osi"].GetBool()
        if (self.calculate_osi):
            wss_stadistics.CalculateOSI(model_part_name,modelpart.process.Info.KratosMultiphysics[step])
            wss_stadistics.CalculateTWSS(model_part_name,modelpart.process.Info.KratosMultiphysics[step])