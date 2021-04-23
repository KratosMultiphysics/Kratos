"""
This module contains an interface to the available response functions
"""
import time as timer
import KratosMultiphysics as Kratos
import KratosMultiphysics.ShapeOptimizationApplication as KSO
from KratosMultiphysics.response_functions.response_function_interface import ResponseFunctionInterface

class TotalVolume(ResponseFunctionInterface):
    def __init__(self, identifier, response_settings, model):
        self.identifier = identifier
        self.response_settings = response_settings

        default_parameters = Kratos.Parameters("""
        {
            "response_type"            : "total_volume",
            "model_part_name"          : "PLEASE_SPECIFY_MODEL_PART_NAME",
            "output_variable_name"     : "PLEASE_PROVIDE_A_3D_VARIABLE_NAME"
        }
        """)

        self.response_settings.ValidateAndAssignDefaults(default_parameters)
        self.model = model

        # this variable is used to store shape sensitivities of total volume w.r.t. shape
        # in a different variable than SHAPE_SENSITIVITY variable. It is a safety measure just in
        # case SHAPE_SENSITIVITY is used in some other response, and not to override them.
        self.variable_name = self.response_settings["output_variable_name"].GetString()
        if (Kratos.KratosGlobals.GetVariableType(self.variable_name) != "Array"):
            raise RuntimeError("Please provide a 3d variable.")

        self.variable = Kratos.KratosGlobals.GetVariable(self.variable_name)

    def Initialize(self):
        self.model_part = self.model[self.response_settings["model_part_name"].GetString()]

    def CalculateValue(self):
        start_time = timer.time()

        self.value = KSO.GeometryUtilities().ComputeVolume(self.model_part)

        Kratos.Logger.PrintInfo(self._GetLabel(), "Time needed for calculating the response value = ",round(timer.time() - start_time,2),"s")

    def CalculateGradient(self):
        start_time = timer.time()

        KSO.GeometryUtilities().ComputeVolumeShapeDerivatives(self.model_part, self.variable)

        Kratos.Logger.PrintInfo(self._GetLabel(), "Time needed for calculating the gradient = ",round(timer.time() - start_time,2),"s")

    def GetValue(self):
        return self.value

    def GetNodalGradient(self, variable):
        if variable != Kratos.SHAPE_SENSITIVITY:
            raise RuntimeError("GetNodalGradient: No gradient for {}!".format(variable.Name))

        gradient = {}
        for node in self.model_part.Nodes:
            gradient[node.Id] = node.GetSolutionStepValue(self.variable)

        return gradient

    @staticmethod
    def _GetLabel():
        return "TotalVolume"


