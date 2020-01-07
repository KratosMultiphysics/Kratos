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
from KratosMultiphysics import Logger
from .response_function import ResponseFunctionBase

class PackagingResponseBase(ResponseFunctionBase):
    """
    A base class for packaging response functions that agglomerate the nodal violations
    into a single response function.
    The agglomeration happens by summing up the square of each nodal violation.
    Nodes that are feasible do NOT contribute to the response value/gradient.
    This is why a prediction of the violation using the gradients is not possible,
    only correction of violations (e.g. from the last step) will happen.

    Derived classes need to implement the calculation of the nodal violations

    Important settings:
    feasible_in_normal_direction : boolean flag that indicates if the normal of bounding instance
        points to the feasible side. True by default
    """

    def __init__(self, identifier, response_settings, model):
        self.identifier = identifier

        response_settings.ValidateAndAssignDefaults(self.GetDefaultSettings())

        self.response_settings = response_settings
        self.model = model

        model_part_name = response_settings["model_part_name"].GetString()
        input_type = response_settings["model_import_settings"]["input_type"].GetString()
        if input_type == "mdpa":
            self.model_part = self.model.CreateModelPart(model_part_name)
            domain_size = response_settings["domain_size"].GetInt()
            if domain_size not in [2, 3]:
                raise Exception("PackagingResponseBase: Invalid 'domain_size': {}".format(domain_size))
            self.model_part.ProcessInfo.SetValue(KM.DOMAIN_SIZE, domain_size)
        elif input_type == "use_input_model_part":
            self.model_part = self.model.GetModelPart(model_part_name)
        else:
            raise Exception("PackagingResponseBase: '{}' model part input type not implemented.".format(input_type))

        self.signed_distances = None
        self.directions = None

        self.value = None

        self.gradient = {}

        self.feasible_in_normal_direction = self.response_settings["feasible_in_normal_direction"].GetBool()
        self.exponent = 2

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
            "feasible_in_normal_direction" : true
        }""")
        return this_defaults

    def Initialize(self):
        if self.response_settings["model_import_settings"]["input_type"].GetString() == "mdpa":
            file_name = self.response_settings["model_import_settings"]["input_filename"].GetString()
            model_part_io = KM.ModelPartIO(file_name)
            model_part_io.ReadModelPart(self.model_part)

    def InitializeSolutionStep(self):
        self.value = None
        self.signed_distances = None
        self.directions = None
        self.gradient = {}

    def CalculateValue(self):
        Logger.PrintInfo("\n> Starting primal analysis for response", self.identifier)

        startTime = timer.time()

        if not self.directions or not self.signed_distances:
            self._CalculateDistances()

        value = 0.0
        for i in range(len(self.signed_distances)):
            value += self._CalculateNodalValue(self.signed_distances[i])

        self.value = value

        Logger.PrintInfo("> Time needed for calculating the response value = ", round(timer.time() - startTime,2), "s")

    def CalculateGradient(self):
        Logger.PrintInfo("\n> Starting gradient calculation for response", self.identifier)

        startTime = timer.time()

        if not self.directions or not self.signed_distances:
            self._CalculateDistances()

        for i, node in enumerate(self.model_part.Nodes):
            gradient = self._CalculateNodalGradient(self.signed_distances[i], self.directions[i*3:i*3+3])
            self.gradient[node.Id] = gradient

        Logger.PrintInfo("> Time needed for calculating gradients = ", round(timer.time() - startTime,2), "s")

    def GetValue(self):
        return self.value

    def GetShapeGradient(self):
        if not self.gradient:
            raise RuntimeError("Gradient was not calculated")
        return self.gradient

    def _CalculateDistances(self):
        raise NotImplementedError("_CalculateDistances needs to be implemented by the derived class!")

    def _CalculateNodalValue(self, signed_distance):
        if not self._HasContribution(signed_distance):
            return 0.0
        return pow(signed_distance, self.exponent)

    def _CalculateNodalGradient(self, signed_distance, direction):
        if not self._HasContribution(signed_distance):
            return [0.0, 0.0, 0.0]

        factor = self.exponent * pow(signed_distance, self.exponent-1)
        gradient = [
            direction[0] * factor,
            direction[1] * factor,
            direction[2] * factor
        ]
        return gradient

    def _HasContribution(self, signed_distance):
        if not self.feasible_in_normal_direction and signed_distance > 0:
            return True
        elif self.feasible_in_normal_direction and signed_distance < 0:
            return True
        return False
