"""
This module contains an interface to the total volume response function
"""
import time as timer
import KratosMultiphysics as Kratos
import KratosMultiphysics.ShapeOptimizationApplication as KSO
try:
    import KratosMultiphysics.StructuralMechanicsApplication as KSM
except ImportError:
    KSM = None
from KratosMultiphysics.response_functions.response_function_interface import ResponseFunctionInterface

class TotalVolume(ResponseFunctionInterface):
    def __init__(self, identifier, response_settings, model):
        self.identifier = identifier
        self.response_settings = response_settings

        default_parameters = Kratos.Parameters("""
        {
            "response_type"            : "total_volume",
            "model_part_name"          : "PLEASE_SPECIFY_MODEL_PART_NAME"
        }
        """)

        self.response_settings.ValidateAndAssignDefaults(default_parameters)
        self.model = model

    def Initialize(self):
        self.model_part = self.model[self.response_settings["model_part_name"].GetString()]

    def CalculateValue(self):
        start_time = timer.time()

        self.value = KSO.GeometryUtilities(self.model_part).ComputeVolume()

        Kratos.Logger.PrintInfo(self._GetLabel(), "Time needed for calculating the response value = ",round(timer.time() - start_time,2),"s")

    def CalculateGradient(self):
        start_time = timer.time()

        KSO.GeometryUtilities(self.model_part).ComputeVolumeShapeDerivatives(Kratos.SHAPE_SENSITIVITY)
        self.gradient = {}
        for node in self.model_part.Nodes:
            self.gradient[node.Id] = node.GetSolutionStepValue(Kratos.SHAPE_SENSITIVITY)

        Kratos.Logger.PrintInfo(self._GetLabel(), "Time needed for calculating the gradient = ",round(timer.time() - start_time,2),"s")

    def GetValue(self):
        return self.value

    def GetNodalGradient(self, variable):
        if variable != Kratos.SHAPE_SENSITIVITY:
            raise RuntimeError("GetNodalGradient: No gradient for {}!".format(variable.Name))
        return self.gradient

    def GetElementalGradient(self, variable):
        if variable != KSM.THICKNESS_SENSITIVITY:
            raise RuntimeError("GetElementalGradient: No gradient for {}!".format(variable.Name))
        gradient = {}
        for condition in self.model_part.Conditions:
            gradient[condition.Id] = 0.0
        return gradient

    @staticmethod
    def _GetLabel():
        return "TotalVolume"
