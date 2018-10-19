import math

import KratosMultiphysics
import KratosMultiphysics.ShapeOptimizationApplication as ShapeOptimizationApplication

class UpperBound(object):
    def __init__(self, bound):
        self.bound = bound

    def CalculateDistance(self, value):
        """negative distance means violation"""
        return self.bound - value

class LowerBound(object):
    def __init__(self, bound):
        self.bound = bound

    def CalculateDistance(self, value):
        """negative distance means violation"""
        return value - self.bound

class ConditionalDampingBase(object):

    def __init__(self, model, settings):
        self.model_part = model[settings["model_part_name"].GetString()]

        self.damping_components = []
        for i in range(3):
            self.damping_components.append(settings["components"][i].GetBool())

    def GetDampingWeights(self, node):
        raise RuntimeError("Called Base class method! Please return [weight_x, weight_y, weight_z]!" )

    def ApplyToNodalSensitivityVariable(self, variable_to_damp):
        self._ApplyToNodalVariable(variable_to_damp)

    def ApplyToNodalUpdateVariable(self, variable_to_damp):
        self._ApplyToNodalVariable(variable_to_damp)

    def _ApplyToNodalVariable(self, variable_to_damp):
        for node in self.model_part.Nodes:

            damping_weights = self.GetDampingWeights(node)

            if damping_weights is None: continue

            vector = node.GetSolutionStepValue(variable_to_damp)
            for i, damping_component in enumerate(self.damping_components):
                if damping_component:
                    vector[i] *= damping_weights[i]

            node.SetSolutionStepValue(variable_to_damp, vector)


