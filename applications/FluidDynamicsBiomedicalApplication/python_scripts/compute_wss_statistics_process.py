# Importing the Kratos Library
from numpy import isin
import KratosMultiphysics

# Import applications
import KratosMultiphysics.FluidDynamicsApplication as KratosFluid
import KratosMultiphysics.FluidDynamicsBiomedicalApplication as KratosBio

def Factory(settings, model):
    if not isinstance(settings, KratosMultiphysics.Parameters):
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
            "calculate_osi": true,
            "compute_normals": true,
            "compute_normals_at_each_step": false
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
        skin_model_part = self.model.GetModelPart(self.settings["model_part_name"].GetString())
        if self.settings["compute_normals"].GetBool():
            KratosMultiphysics.NormalCalculationUtils().CalculateOnSimplexNonHistorical(
                skin_model_part,
                skin_model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE],
                KratosFluid.WALL_NORMAL)

        # Initialize the WSS variables
        KratosBio.WssStatisticsUtilities.InitializeWSSVariables(skin_model_part)

    def ExecuteBeforeOutputStep(self):
        skin_model_part = self.model.GetModelPart(self.settings["model_part_name"].GetString())

        if self.settings["compute_normals_at_each_step"].GetBool():
            KratosMultiphysics.NormalCalculationUtils().CalculateOnSimplexNonHistorical(
                skin_model_part,
                skin_model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE],
                KratosFluid.WALL_NORMAL)

        if self.settings["calculate_wss"].GetBool():
            is_normal_historical = False
            KratosBio.WssStatisticsUtilities.CalculateWSS(skin_model_part, KratosFluid.WALL_NORMAL, is_normal_historical)

        if self.settings["calculate_osi"].GetBool():
            KratosBio.WssStatisticsUtilities.CalculateOSI(skin_model_part)
