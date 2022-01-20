# importing the Kratos Library
import KratosMultiphysics as KM
from KratosMultiphysics import Parameters, Logger
from KratosMultiphysics.response_functions.response_function_interface import ResponseFunctionInterface
import KratosMultiphysics.ShapeOptimizationApplication as KSO

import time as timer
import math


def _AddConditionsFromParent(parent, child):
    node_ids = set([node.Id for node in child.Nodes])
    conditions = []
    for condition in parent.Conditions:
        all_nodes_found = True
        for node in condition.GetNodes():
            if node.Id not in node_ids:
                all_nodes_found = False
                break
        if all_nodes_found:
            conditions.append(condition.Id)
    child.AddConditions(conditions)

# ==============================================================================
class ConstantResponseFunction(ResponseFunctionInterface):
    """
    Constant value response function.

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
                raise Exception("ConstantResponseFunction: Invalid 'domain_size': {}".format(domain_size))
            self.model_part.ProcessInfo.SetValue(KM.DOMAIN_SIZE, domain_size)
            self.model_part_needs_to_be_imported = True
        elif input_type == "use_input_model_part":
            self.model_part = None  # will be retrieved in Initialize()
        else:
            raise Exception("Other model part input options are not yet implemented.")

        self.response_function_utility = None  # will be created in Initialize()

        self.model.GetModelPart(self._model_part_name.split(".")[0]).AddNodalSolutionStepVariable(KM.SHAPE_SENSITIVITY)

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
            "values": [1.0, 1.0, 1.0]
        }""")
        return this_defaults

    def Initialize(self):
        if self.response_settings["model_import_settings"]["input_type"].GetString() == "mdpa":
            file_name = self.response_settings["model_import_settings"]["input_filename"].GetString()
            model_part_io = KM.ModelPartIO(file_name)
            model_part_io.ReadModelPart(self.model_part)
        else:
            self.model_part = self.model.GetModelPart(self._model_part_name)

        only = self.response_settings["only"].GetString()
        if only != "":
            only_part = self.model.GetModelPart(only)
            if only_part.NumberOfConditions() == 0:
                _AddConditionsFromParent(self.model_part, only_part)
                Logger.PrintWarning("ConstantResponse", "Automatically added {} conditions to model_part '{}'.".format(only_part.NumberOfConditions(), only_part.Name))
        else:
            only_part = self.model_part

        if only_part.NumberOfConditions() == 0:
            raise RuntimeError("The model_part '{}' does not have any surface conditions!".format(only_part.Name))


    def UpdateDesign(self, updated_model_part, variable):
        self.value = None

    def CalculateValue(self):
        Logger.PrintInfo("ConstantResponse", "Starting calculation of response value:", self.identifier)

        startTime = timer.time()
        self.value = 1
        Logger.PrintInfo("ConstantResponse", "Time needed for calculating the response value = ",round(timer.time() - startTime,2),"s")

    def CalculateGradient(self):
        Logger.PrintInfo("ConstantResponse", "Starting gradient calculation for response", self.identifier)

        for node in self.model_part.Nodes:
            area_normal = node.GetSolutionStepValue(KM.NORMAL, 0)
            area = math.sqrt(area_normal[0]**2 + area_normal[1]**2 + area_normal[2]**2)
            val = 1 * area             
            node.SetSolutionStepValue(KM.SHAPE_SENSITIVITY, [0, 0, val])

        startTime = timer.time()
        Logger.PrintInfo("ConstantResponse", "Time needed for calculating gradients",round(timer.time() - startTime,2),"s")

    def GetValue(self):
        return self.value

    def GetNodalGradient(self, variable):
        if variable != KM.SHAPE_SENSITIVITY:
            raise RuntimeError("GetNodalGradient: No gradient for {}!".format(variable.Name))
        gradient = {}
        for node in self.model_part.Nodes:
            gradient[node.Id] = node.GetSolutionStepValue(variable)
        return gradient
