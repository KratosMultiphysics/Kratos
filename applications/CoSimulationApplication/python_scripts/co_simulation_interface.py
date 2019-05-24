import numpy as np

import KratosMultiphysics.CoSimulationApplication.co_simulation_tools as cs_tools
cs_data_structure = cs_tools.cs_data_structure


# Class CoSimulationInterface: Holds the different ModelParts of the interface.
class CoSimulationInterface(object):
    def __init__(self, model_parts):
        super().__init__()

        self.model_parts = model_parts

    def GetPythonList(self):
        data = []
        step = 0
        for model_part in self.model_parts:
            for node in model_part.Nodes:
                for variable in list(node.variables[step].keys()):
                    value = node.GetSolutionStepValue(variable, step)
                    data.append(value)
        return data

    def GetNumpyArray(self):
        return np.array(self.GetPythonList())

    def SetPythonList(self, data):
        index = 0
        step = 0
        for model_part in self.model_parts:
            for node in model_part.Nodes:
                for variable in list(node.variables[step].keys()):
                    value = data[index]
                    index += 1
                    node.SetSolutionStepValue(variable, step, value)

    def SetNumpyArray(self, data):
        self.SetPythonList(data.tolist())

    def __add__(self, other):
        # To do: use zip for looping over ModelParts of two interfaces simultaneously.
        print(type(self), type(other))
        return other

    def __sub__(self, other):
        print(type(self), type(other))
        return other

    def __mul__(self, other):
        print(type(self), type(other))
        return self
