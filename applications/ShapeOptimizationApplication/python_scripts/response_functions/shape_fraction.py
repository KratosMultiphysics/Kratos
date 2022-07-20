# importing the Kratos Library
import KratosMultiphysics as KM
from KratosMultiphysics import Parameters, Logger
from KratosMultiphysics.response_functions.response_function_interface import ResponseFunctionInterface
import KratosMultiphysics.ShapeOptimizationApplication as KSO
from KratosMultiphysics.ShapeOptimizationApplication.utilities import custom_math as cm
from KratosMultiphysics.ShapeOptimizationApplication.utilities.custom_variable_utilities import ReadNodalVariableToList

import time as timer

# ==============================================================================
class ShapeFractionResponseFunction(ResponseFunctionInterface):
    """
    Shape Fraction response function.
    It sums up all the nodes which have a shape change.
    Nodes that are feasible do NOT contribute to the response value/gradient.
    This is why a prediction of the violation using the gradients is not possible,
    only correction of violations (e.g. from the last step) will happen.

    Attributes
    ----------
    model_part : Model part object of the response function
    response_function_utility: Cpp utilities object doing the actual computation of response value and gradient.
    """

    def __init__(self, identifier, response_settings, model):
        self.identifier = identifier

        response_settings.ValidateAndAssignDefaults(self.GetDefaultParameters())

        self.response_settings = response_settings
        self.model = model
        self.model_part_needs_to_be_imported = False

        self.value = None

        self._model_part_name = response_settings["model_part_name"].GetString()
        input_type = response_settings["model_import_settings"]["input_type"].GetString()
        if input_type == "mdpa":
            self.model_part = self.model.CreateModelPart(self._model_part_name, 2)
            domain_size = response_settings["domain_size"].GetInt()
            if domain_size not in [3]:
                raise Exception("ShapeFractionResponseFunction: Invalid 'domain_size': {}".format(domain_size))
            self.model_part.ProcessInfo.SetValue(KM.DOMAIN_SIZE, domain_size)
            self.model_part_needs_to_be_imported = True
        elif input_type == "use_input_model_part":
            self.model_part = None  # will be retrieved in Initialize()
        else:
            raise Exception("Other model part input options are not yet implemented.")

        self.tolerance = response_settings["tolerance"].GetDouble()
        self.max_shape_fraction = response_settings["max_shape_fraction"].GetDouble()
        #self.model.GetModelPart(self._model_part_name.split(".")[0]).AddNodalSolutionStepVariable(KM.SHAPE_SENSITIVITY)

    @classmethod
    def GetDefaultParameters(cls):
        this_defaults = KM.Parameters("""{
            "response_type"         : "UNKNOWN_TYPE",
            "model_part_name"       : "UNKNOWN_NAME",
            "only"                  : "",
            "domain_size"           : 3,
            "model_import_settings" : {
                "input_type"        : "use_input_model_part",
                "input_filename"    : "UNKNOWN_NAME"
            },
            "max_shape_fraction": 0.5,
            "tolerance": 1e-6
        }""")
        return this_defaults

    def Initialize(self):
        if self.response_settings["model_import_settings"]["input_type"].GetString() == "mdpa":
            file_name = self.response_settings["model_import_settings"]["input_filename"].GetString()
            model_part_io = KM.ModelPartIO(file_name)
            model_part_io.ReadModelPart(self.model_part)
        else:
            self.model_part = self.model.GetModelPart(self._model_part_name)

    def InitializeSolutionStep(self):
        self.value = None
        self.gradient = {}

    def CalculateValue(self):
        Logger.PrintInfo("ShapeFractionResponse", "Starting calculation of response value:", self.identifier)

        startTime = timer.time()
        self.value = 0.0
        total_nodes = len(self.model_part.Nodes)
        g = 0.0
        for node in self.model_part.Nodes:
            g_i = self._CalculateNodalValue(node)
            g += g_i
        self.value = g / total_nodes - self.max_shape_fraction

        Logger.PrintInfo("ShapeFractionResponse", "Time needed for calculating the response value = ",round(timer.time() - startTime,2),"s")

    def CalculateGradient(self):
        Logger.PrintInfo("ShapeFractionResponse", "Starting gradient calculation for response", self.identifier)

        startTime = timer.time()
        # quantile = self._GetQuantile()
        quantile = 1.0
        for node in self.model_part.Nodes:
            gradient_i = self._CalculateNodalGradient(node, quantile)
            self.gradient[node.Id] = gradient_i

        Logger.PrintInfo("ShapeFractionResponse", "Time needed for calculating gradients",round(timer.time() - startTime,2),"s")

    def GetValue(self):
        return self.value

    def GetNodalGradient(self, variable):
        if not self.gradient:
            raise RuntimeError("Gradient was not calculated")
        if variable != KM.SHAPE_SENSITIVITY:
            raise RuntimeError("GetNodalGradient: No gradient for {}!".format(variable.Name))
        return self.gradient

    def _CalculateNodalValue(self, node):
        shape_change = node.GetSolutionStepValue(KSO.SHAPE_CHANGE)
        norm = cm.Norm2(shape_change)

        if norm > self.tolerance:
            return 1.0
        else:
            return 0.0

    def _CalculateNodalGradient(self, node, quantile):
        shape_change = node.GetSolutionStepValue(KSO.SHAPE_CHANGE)
        norm = cm.Norm2(shape_change)

        if norm > self.tolerance:
            return [
                (shape_change[0] / norm)*(quantile / norm)**1,
                (shape_change[1] / norm)*(quantile / norm)**1,
                (shape_change[2] / norm)*(quantile / norm)**1
            ]
        else:
            return [0.0, 0.0, 0.0]

    def _GetQuantile(self):

        nodal_variable = KM.KratosGlobals.GetVariable("SHAPE_CHANGE")
        shape_change = ReadNodalVariableToList(self.model_part, nodal_variable)

        shape_change_norm = []
        Logger.PrintInfo("print(len(shape_change)/3):{}".format(len(shape_change)/3))
        for i in range(int(len(shape_change)/3)):
            shape_change_norm.append(cm.Norm2(shape_change[3*i:3*i+2]))

        shape_change_norm.sort()

        index = round(len(shape_change_norm)*(1-self.max_shape_fraction))

        quantile = shape_change_norm[index]

        Logger.PrintInfo("ShapeChange Quantile {}: {}".format((1-self.max_shape_fraction), quantile))

        return quantile
