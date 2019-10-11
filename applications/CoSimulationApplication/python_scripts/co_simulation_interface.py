import numpy as np
import copy

import KratosMultiphysics.CoSimulationApplication.co_simulation_tools as cs_tools
cs_data_structure = cs_tools.cs_data_structure
import KratosMultiphysics as KM


# Class CoSimulationInterface: Holds the different ModelParts of the interface.
class CoSimulationInterface(object):
    """
    When copying a CosimulationInterface object, its proper
    deepcopy method must be used (self.deepcopy()). If they are
    copied using copy.deepcopy(), the global Variables will be
    referenced wrong.
    """

    def __init__(self, model, parameters):
        """
        The input 'parameters' is a Parameters object: the keys
        contain the ModelPart names, the values contain the
        names of the Variables. The latter are given either as
        a single string, or as a list of one or more strings.

        Two examples of input Parameters:
            {
            "interface_a":
                {
                    "mp_1": ["PRESSURE", "FORCE"],
                    "mp_2": "DENSITY"
                }
            }
            {
            "interface_b":
                {
                    "mp_1": ["TEMPERATURE"],
                    "mp_3": "DISPLACEMENT"
                }
            }
        """

        super().__init__()

        self.model = model
        self.model_parts_variables = list(parameters.items())

    def GetPythonList(self):
        data = []
        step = 0
        for model_part_name, variable_names in self.model_parts_variables:
            model_part = self.model.GetModelPart(model_part_name)
            for variable_name in variable_names.list():
                variable = vars(KM)[variable_name.GetString()]
                for node in model_part.Nodes:
                    value = node.GetSolutionStepValue(variable, step)
                    if type(value) is list:
                        data += value
                    else:
                        data.append(value)
        return data

    def GetNumpyArray(self):
        return np.array(self.GetPythonList())

    def SetPythonList(self, data):
        index = 0
        step = 0
        for model_part_name, variable_names in self.model_parts_variables:
            model_part = self.model.GetModelPart(model_part_name)
            for variable_name in variable_names.list():
                variable = vars(KM)[variable_name.GetString()]
                for node in model_part.Nodes:
                    if variable.Type() is "Double":
                        value = data[index]
                        index += 1
                    elif variable.Type() is "Array":
                        value = data[index:index + 3]
                        index += 3
                    else:
                        raise NotImplementedError('Only "Double" and "Array" Variables implemented.')
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
        # *** --> might be problematic, because method calls itself...
        # *** in that case, I should do initial deepcopy manually or sth...

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
        if isinstance(other, CoSimulationInterface):
            result.SetNumpyArray(self.GetNumpyArray() + other.GetNumpyArray())
        elif isinstance(other, (int, float, np.integer, np.floating)):
            result.SetNumpyArray(self.GetNumpyArray() + other)
        else:
            return NotImplemented
        return result

    def __radd__(self, other):
        return self.__add__(other)

    def __sub__(self, other):
        result = self.deepcopy()
        if isinstance(other, CoSimulationInterface):
            result.SetNumpyArray(self.GetNumpyArray() - other.GetNumpyArray())
        elif isinstance(other, (int, float, np.integer, np.floating)):
            result.SetNumpyArray(self.GetNumpyArray() - other)
        else:
            return NotImplemented
        return result

    def __rsub__(self, other):
        return self.__sub__(other).__mul__(-1)

    def __mul__(self, other):
        result = self.deepcopy()
        if isinstance(other, CoSimulationInterface):
            result.SetNumpyArray(self.GetNumpyArray() * other.GetNumpyArray())
        elif isinstance(other, (int, float, np.integer, np.floating)):
            result.SetNumpyArray(self.GetNumpyArray() * other)
        else:
            return NotImplemented
        return result

    def __rmul__(self, other):
        return self.__mul__(other)

    def __truediv__(self, other):
        result = self.deepcopy()
        if isinstance(other, (int, float, np.integer, np.floating)):
            result.SetNumpyArray(self.GetNumpyArray() / other)
        else:
            return NotImplemented
        return result

    def __pow__(self, other):
        result = self.deepcopy()
        if isinstance(other, (int, float, np.integer, np.floating)):
            result.SetNumpyArray(self.GetNumpyArray() ** other)
        else:
            return NotImplemented
        return result
