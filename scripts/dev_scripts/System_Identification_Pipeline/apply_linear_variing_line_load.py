from math import sin, pi

import KratosMultiphysics as Kratos
import KratosMultiphysics.StructuralMechanicsApplication as KratosSM


def Factory(parameters, model):
    return ApplyLinearVariingLineLoad(model, parameters)


class ApplyLinearVariingLineLoad(Kratos.Process):
    def __init__(self, model, parameters):
        super().__init__()
        parameters = parameters["Parameters"]
        self.model_part = model[parameters["model_part_name"].GetString()]
        self.strength = parameters["modulus"].GetDouble()

    def ExecuteInitializeSolutionStep(self):

        # Get hight of model part
        y_max = -10e64
        y_min = +10e64
        for node in self.model_part.Nodes:
            if node.Y > y_max:
                y_max = node.Y
            if node.Y < y_min:
                y_min = node.Y

        def line_load_strength(y, y_min, y_max, strength):
            return strength * (2 * (y-((y_min+y_max)/2))) / (y_max - y_min)

        conditions = self.model_part.GetConditions()

        for condition in conditions:

            # Get nodes of condition
            nodes = condition.GetNodes()
            upper_node = nodes[0]
            lower_node = nodes[1]
            if lower_node.Y > upper_node.Y:
                upper_node = nodes[1]
                lower_node = nodes[0]

            # find pressure for line load segment
            strength_upper_node = line_load_strength(upper_node.Y, y_min=y_min, y_max=y_max, strength=self.strength)
            strength_lower_node = line_load_strength(lower_node.Y, y_min=y_min, y_max=y_max, strength=self.strength)

            # The nodes at the sides (max/min height) only get half of the force
            if upper_node.Y == y_max:
                strength_upper_node /= 2.0

            if lower_node.Y == y_min:
                strength_lower_node /= 2.0

            # Apply force to nodes - requires that there is a line load condition in the mdpa at the nodes that are loaded
            upper_node.SetSolutionStepValue(KratosSM.LINE_LOAD, Kratos.Array3([strength_upper_node, 0.0, 0.0]))
            lower_node.SetSolutionStepValue(KratosSM.LINE_LOAD, Kratos.Array3([strength_lower_node, 0.0, 0.0]))
