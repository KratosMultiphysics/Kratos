import math

import KratosMultiphysics
import KratosMultiphysics.ShapeOptimizationApplication as ShapeOptimizationApplication
import conditional_damping_base

def Factory(model, settings):
    if not isinstance(settings, KratosMultiphysics.Parameters):
        raise Exception("Expected input shall be a Parameters object, encapsulating a json string")
    return CoordinateBoundDamping(model, settings["Parameters"])


class CoordinateBoundDamping(conditional_damping_base.ConditionalDampingBase):

    def __init__(self, model, settings):

        default_settings = KratosMultiphysics.Parameters(
            """
            {
                "model_part_name" : "SPECIFY_MODEL_PART_NAME",
                "type"            : "",
                "condition"       : "",
                "value"           : -1e30,
                "transition_size" : 0.0,
                "components"      : [true, true, true]
            }
            """
        )

        settings.ValidateAndAssignDefaults(default_settings)

        super(CoordinateBoundDamping, self).__init__(model, settings)

        self.settings = settings

        self.transition_size = self.settings["transition_size"].GetDouble()
        if self.transition_size < 0.0:
            raise RuntimeError("Specify 'transition_size' >= 0.0!")

        _value = self.settings["value"].GetDouble()
        if _value == 1e30:
            raise RuntimeError("Compulsory settings 'value' is missing!")

        _type = self.settings["type"].GetString()
        _available_types = ["X", "Y", "Z"]

        if _value == 1e30:
            raise RuntimeError("Compulsory settings 'value' is missing!")

        _condition = self.settings["condition"].GetString()
        _available_conditions = ["<", ">"]
        if _condition not in _available_conditions:
            raise RuntimeError("Settings 'condition' is missing or invalid! Possible settings are: " + str(_available_conditions))

        elif _type == "X":
            self.GetValue = lambda node : node.X0
        elif _type == "Y":
            self.GetValue = lambda node : node.Y0
        elif _type == "Z":
            self.GetValue = lambda node : node.Z0
        else:
            raise RuntimeError("Settings 'type' is missing or invalid! Possible settings are: " + str(_available_types))

        if _condition == "<":
            self.bound = conditional_damping_base.UpperBound(_value)
        else:
            self.bound = conditional_damping_base.LowerBound(_value)

    def GetDampingWeights(self, node):
        scalar_value = self.GetValue(node)
        distance = self.bound.CalculateDistance(scalar_value)
        if distance < 0.0:
            damping_weight = 0.0
        elif distance < self.transition_size:
            damping_weight = distance / self.transition_size
        else:
            return None

        return [damping_weight, damping_weight, damping_weight]


