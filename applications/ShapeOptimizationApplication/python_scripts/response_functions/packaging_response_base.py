from __future__ import print_function, absolute_import, division

import time as timer

import KratosMultiphysics
from KratosMultiphysics import Logger

from .response_function import ResponseFunctionBase

import numpy as np

class ValueAndGradient(object):

    def __init__(self, infeasible_side):
        self.infeasible_side = infeasible_side

    def GetValue(self, distance, normal):
        if self.infeasible_side and distance > 0:
            return self._GetValue(distance, normal)
        elif not self.infeasible_side and distance < 0:
            return self._GetValue(distance, normal)
        return 0.0

    def GetGradient(self, distance, normal):
        if self.infeasible_side and distance > 0:
            return self._GetGradient(distance, normal)
        elif not self.infeasible_side and distance < 0:
            return self._GetGradient(distance, normal)
        return [0.0, 0.0, 0.0]

    def _GetValue(self, distance, normal):
        raise NotImplementedError()

    def _GetGradient(self, distance, normal):
        raise NotImplementedError()

class ScalarDistancePow(ValueAndGradient):

    exponent = 1

    def _GetValue(self, distance, normal):
        return pow(distance, self.exponent)

    def _GetGradient(self, distance, normal):
        factor = self.exponent * pow(distance, self.exponent-1)
        gradient = [
            normal[0] * factor,
            normal[1] * factor,
            normal[2] * factor
        ]
        return gradient

class ComponentDistance(ValueAndGradient):
    '''Is the same as ScalarDistancePow with exponent 1'''

    def _GetValue(self, distance, normal):
        return distance

    def _GetGradient(self, distance, normal):
        gradient = [
            normal[0],
            normal[1],
            normal[2]
        ]
        return gradient

class Exponential(object):

    def _GetValue(self, distance, normal):
        return np.exp(distance)

    def _GetGradient(self, distance, normal):
        factor = np.exp(distance)
        gradient = [
            normal[0] * factor,
            normal[1] * factor,
            normal[2] * factor
        ]
        return gradient

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
                raise Exception("MassResponseFunction: Invalid 'domain_size': {}".format(domain_size))
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

        method_class = globals()[self.response_settings["method"].GetString()]
        self.method = method_class(self.response_settings["infeasible_side"].GetBool())

    def Initialize(self):
        if self.model_part_needs_to_be_imported:
            # import model part
            model_part_io = KratosMultiphysics.ModelPartIO(self.response_settings["model_import_settings"]["input_filename"].GetString())
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

            value += self.method.GetValue(self.signed_distances[i], _direction)

        self.value = value

        print(value)
        print("#"*80)
        Logger.PrintInfo("> Time needed for calculating the response value = ",round(timer.time() - startTime,2),"s")

    def CalculateGradient(self):
        Logger.PrintInfo("\n> Starting gradient calculation for response", self.identifier)

        startTime = timer.time()

        if not self.directions or not self.signed_distances:
            self._CalculateProjectedDistances()

        i = 0
        for node in self.model_part.Nodes:

            _direction = np.array(self.directions[i*3:i*3+3])

            gradient = self.method.GetGradient(self.signed_distances[i], _direction)

            self.gradient[node.Id] = gradient

            i+=1

        Logger.PrintInfo("> Time needed for calculating gradients",round(timer.time() - startTime,2),"s")

    def GetValue(self):
        return self.value

    def GetShapeGradient(self):
        if not self.gradient:
            raise RuntimeError("Gradient was not calculated")
        return self.gradient

    def _CalculateProjectedDistances(self):
        raise NotImplementedError("_CalculateProjectedDistances needs to be implemented by the derived class!")

