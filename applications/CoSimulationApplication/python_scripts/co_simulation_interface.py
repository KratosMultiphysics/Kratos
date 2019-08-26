import numpy as np
import copy

import KratosMultiphysics.CoSimulationApplication.co_simulation_tools as cs_tools
cs_data_structure = cs_tools.cs_data_structure


# Class CoSimulationInterface: Holds the different ModelParts of the interface.
class CoSimulationInterface(object):
    def __init__(self, model, parameters):
        super().__init__()

        self.model = model
        self.model_parts_variables = list(parameters.items())

    def GetPythonList(self):
        data = []
        step = 0
        for model_part_name, variable in self.model_parts_variables:
            model_part = self.model.GetModelPart(model_part_name)
            for node in model_part.Nodes:
                value = node.GetSolutionStepValue(variable, step)
                data.append(value)
        return data

    def GetNumpyArray(self):
        return np.array(self.GetPythonList())

    def SetPythonList(self, data):
        index = 0
        step = 0
        for model_part_name, variable in self.model_parts_variables:
            model_part = self.model.GetModelPart(model_part_name)
            for node in model_part.Nodes:
                value = data[index]
                index += 1
                node.SetSolutionStepValue(variable, step, value)

    def SetNumpyArray(self, data):
        self.SetPythonList(data.tolist())

    def __add__(self, other):
        result = copy.deepcopy(self)
        result.SetNumpyArray(self.GetNumpyArray() + other.GetNumpyArray())
        return result

    def __sub__(self, other):
        result = copy.deepcopy(self)
        result.SetNumpyArray(self.GetNumpyArray() - other.GetNumpyArray())
        return result

    def __mul__(self, other):
        if type(other) == float:
            result = copy.deepcopy(self)
            result.SetNumpyArray(self.GetNumpyArray() * other)
            return result
        else:
            Exception("Not implemented.")
