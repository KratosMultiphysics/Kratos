import numpy as np
import copy

import KratosMultiphysics.CoSimulationApplication.co_simulation_tools as cs_tools
cs_data_structure = cs_tools.cs_data_structure
import KratosMultiphysics as KM


# Class CoSimulationInterface: Holds the different ModelParts of the interface.
class CoSimulationInterface(object):
    def __init__(self, model, parameters):
        super().__init__()

        self.model = model
        self.model_parts_variables = list(parameters.items())

    def GetPythonList(self):
        data = []
        step = 0
        for model_part_name, variable_name in self.model_parts_variables:
            model_part = self.model.GetModelPart(model_part_name)
            variable = vars(KM)[variable_name.GetString()]
            for node in model_part.Nodes:
                value = node.GetSolutionStepValue(variable, step)
                data.append(value)
        return data

    def GetNumpyArray(self):
        return np.array(self.GetPythonList())

    def SetPythonList(self, data):
        index = 0
        step = 0
        for model_part_name, variable_name in self.model_parts_variables:
            model_part = self.model.GetModelPart(model_part_name)
            variable = vars(KM)[variable_name.GetString()]
            for node in model_part.Nodes:
                value = data[index]
                index += 1
                node.SetSolutionStepValue(variable, step, value)

    def SetNumpyArray(self, data):
        self.SetPythonList(data.tolist())

    def deepcopy(self):
        """
        WARNING: current implementation of this method works only for
        simple Model, not if it contains nested ModelParts!

        It is not allowed to use the copy.deepcopy() for copying
        CoSimulationInterface objects.
        The reason is that the __hist_variables (AREA, PRESSURE, ...)
        are also copied. However, these are global variables, so
        when they're copied, their address changes which gives problems.
        """

        # *** TODO: change method to __deepcopy__(self, _) and test it

        cp = copy.deepcopy(self)
        cp.model_parts_variables = self.model_parts_variables
        for key_mp in cp.model.__dict__.keys():
            cp.model.GetModelPart(key_mp)._ModelPart__hist_variables = \
                self.model.GetModelPart(key_mp)._ModelPart__hist_variables
            # print([id(_) for _ in model.GetModelPart(key)._ModelPart__hist_variables])
            # print([id(_) for _ in self.model.GetModelPart(key)._ModelPart__hist_variables])
            for key_n in cp.model.GetModelPart(key_mp)._ModelPart__nodes.keys():
                cp.model.GetModelPart(key_mp)._ModelPart__nodes[key_n]._Node__hist_variables = \
                    self.model.GetModelPart(key_mp)._ModelPart__nodes[key_n]._Node__hist_variables
                # print([id(_) for _ in cp.model.GetModelPart(key_mp)._ModelPart__nodes[key_n]._Node__hist_variables])
                # print([id(_) for _ in self.model.GetModelPart(key_mp)._ModelPart__nodes[key_n]._Node__hist_variables])
                node_n = cp.model.GetModelPart(key_mp)._ModelPart__nodes[key_n]
                node_o = self.model.GetModelPart(key_mp)._ModelPart__nodes[key_n]
                new_dict = {}
                for item in node_o._Node__solution_steps_nodal_data.items():
                    new_dict[item[0]] = copy.deepcopy(item[1])
                node_n._Node__solution_steps_nodal_data = new_dict
                # print([id(_) for _ in node_n._Node__solution_steps_nodal_data.keys()])
                # print([id(_) for _ in node_o._Node__solution_steps_nodal_data.keys()])
        return cp

    def __add__(self, other):
        result = self.deepcopy()
        result.SetNumpyArray(self.GetNumpyArray() + other.GetNumpyArray())
        return result

    def __sub__(self, other):
        result = self.deepcopy()
        result.SetNumpyArray(self.GetNumpyArray() - other.GetNumpyArray())
        return result

    def __mul__(self, other):
        if type(other) == float:
            result = self.deepcopy()
            result.SetNumpyArray(self.GetNumpyArray() * other)
            return result
        else:
            Exception("Not implemented.")
