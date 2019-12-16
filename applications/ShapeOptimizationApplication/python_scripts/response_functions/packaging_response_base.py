from __future__ import print_function, absolute_import, division

import time as timer

import KratosMultiphysics
from KratosMultiphysics import Logger

from .response_function import ResponseFunctionBase


class PackagingResponseBase(ResponseFunctionBase):
    """

    Attributes
    ----------
    model_part : Model part object of the response function
    # TODO
    """

    def __init__(self, identifier, response_settings, model):
        self.identifier = identifier

        self.response_settings = response_settings
        self.model = model
        self.model_part_needs_to_be_imported = False

        model_part_name = response_settings["model_part_name"].GetString()
        input_type = response_settings["model_import_settings"]["input_type"].GetString()
        if input_type == "mdpa":
            self.model_part = self.model.CreateModelPart(model_part_name, 2)
            domain_size = response_settings["domain_size"].GetInt()
            if domain_size not in [2, 3]:
                raise Exception("PackagingResponseBase: Invalid 'domain_size': {}".format(domain_size))
            self.model_part.ProcessInfo.SetValue(KratosMultiphysics.DOMAIN_SIZE, domain_size)
            self.model_part_needs_to_be_imported = True
        elif input_type == "use_input_model_part":
            self.model_part = self.model.GetModelPart(model_part_name)
        else:
            raise Exception("Other model part input options are not yet implemented.")

        self.signed_distances = None
        self.directions = None

        self.value = None

        self.gradient = {}

        self.infeasible_side = self.response_settings["infeasible_side"].GetBool()
        self.exponent = 2

    def Initialize(self):
        if self.model_part_needs_to_be_imported:
            file_name = self.response_settings["model_import_settings"]["input_filename"].GetString()
            model_part_io = KratosMultiphysics.ModelPartIO(file_name)
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
            self._CalculateProjectedDistances()

        value = 0.0
        for i in range(len(self.signed_distances)):

            _direction = np.array(self.directions[i*3:i*3+3])

            value += self._CalculateValueForNode(self.signed_distances[i], _direction)

        self.value = value

        Logger.PrintInfo("> Time needed for calculating the response value = ", round(timer.time() - startTime,2), "s")

    def CalculateGradient(self):
        Logger.PrintInfo("\n> Starting gradient calculation for response", self.identifier)

        startTime = timer.time()

        if not self.directions or not self.signed_distances:
            self._CalculateProjectedDistances()

        i = 0
        for node in self.model_part.Nodes:

            _direction = np.array(self.directions[i*3:i*3+3])

            gradient = self._CalculateGradientForNode(self.signed_distances[i], _direction)

            self.gradient[node.Id] = gradient

            i+=1

        Logger.PrintInfo("> Time needed for calculating gradients = ", round(timer.time() - startTime,2), "s")

    def GetValue(self):
        return self.value

    def GetShapeGradient(self):
        if not self.gradient:
            raise RuntimeError("Gradient was not calculated")
        return self.gradient

    def _CalculateProjectedDistances(self):
        raise NotImplementedError("_CalculateProjectedDistances needs to be implemented by the derived class!")

    def _CalculateValueForNode(self, signed_distance, direction):
        if not self._HasContribution(signed_distance):
            return 0.0
        return pow(signed_distance, self.exponent)

    def _CalculateGradientForNode(self, signed_distance, direction):
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
        if self.infeasible_side and signed_distance > 0:
            return True
        elif not self.infeasible_side and signed_distance < 0:
            return True
        return False


