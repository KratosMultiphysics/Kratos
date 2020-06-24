# ==============================================================================
#  KratosShapeOptimizationApplication
#
#  License:         BSD License
#                   license: ShapeOptimizationApplication/license.txt
#
#  Main authors:    Geiser Armin, https://github.com/armingeiser
#
# ==============================================================================

import time as timer
import KratosMultiphysics as KM
import KratosMultiphysics.ShapeOptimizationApplication as KSO
from KratosMultiphysics import Logger
from KratosMultiphysics.response_functions.response_function_interface import ResponseFunctionInterface


class SurfaceArea(ResponseFunctionInterface):
    """
    This simple geometric response function calculates surface area based on surface conditions"""

    def __init__(self, identifier, response_settings, model):
        self.identifier = identifier

        response_settings.ValidateAndAssignDefaults(self.GetDefaultSettings())

        self.model = model

        self.model_part_needs_to_be_imported = False
        model_part_name = response_settings["model_part_name"].GetString()
        input_type = response_settings["model_import_settings"]["input_type"].GetString()
        if input_type == "mdpa":
            self.model_part = self.model.CreateModelPart(model_part_name)
            domain_size = response_settings["domain_size"].GetInt()
            if domain_size not in [2, 3]:
                raise Exception("SurfaceArea: Invalid 'domain_size': {}".format(domain_size))
            self.model_part.ProcessInfo.SetValue(KM.DOMAIN_SIZE, domain_size)
            self.model_part_needs_to_be_imported = True
        elif input_type == "use_input_model_part":
            self.model_part = self.model.GetModelPart(model_part_name)
        else:
            raise Exception("SurfaceArea: '{}' model part input type not implemented.".format(input_type))

        self.value = None

        self.gradient = {}

        self.response_function_utility = KSO.SurfaceAreaResponseFunctionUtility(self.model_part, response_settings)

    @classmethod
    def GetDefaultSettings(cls):
        this_defaults = KM.Parameters("""{
            "response_type"         : "UNKNOWN_TYPE",
            "model_part_name"       : "UNKNOWN_NAME",
            "domain_size"           : 3,
            "model_import_settings" : {
                "input_type"        : "use_input_model_part",
                "input_filename"    : "UNKNOWN_NAME"
            },
            "perturbation_size" : 1e-6
        }""")
        return this_defaults

    def Initialize(self):
        if self.model_part_needs_to_be_imported:
            model_part_io = KM.ModelPartIO(self.response_settings["model_import_settings"]["input_filename"].GetString())
            model_part_io.ReadModelPart(self.model_part)

    def CalculateValue(self):
        Logger.PrintInfo("SurfaceArea", "Starting primal analysis for response", self.identifier)

        startTime = timer.time()
        self.value = self.response_function_utility.CalculateValue()
        Logger.PrintInfo("SurfaceArea", "Time needed for calculating the response value = ",round(timer.time() - startTime,2),"s")

    def CalculateGradient(self):
        Logger.PrintInfo("SurfaceArea", "Starting gradient calculation for response", self.identifier)

        startTime = timer.time()
        self.response_function_utility.CalculateGradient()
        Logger.PrintInfo("SurfaceArea", "Time needed for calculating gradients",round(timer.time() - startTime,2),"s")

    def GetValue(self):
        return self.value

    def GetNodalGradient(self, variable):
        if variable != KM.SHAPE_SENSITIVITY:
            raise RuntimeError("GetNodalGradient: No gradient for {}!".format(variable.Name))
        gradient = {}
        for node in self.model_part.Nodes:
            gradient[node.Id] = node.GetSolutionStepValue(variable)
        return gradient
