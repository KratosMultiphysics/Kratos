import KratosMultiphysics.CoSimulationApplication.co_simulation_tools as cs_tools
cs_data_structure = cs_tools.cs_data_structure


# Class CoSimulationInterface: Holds the different ModelParts of the interface.
class CoSimulationInterface(object):
    def __init__(self, model_parts):
        super().__init__()

        self.model_parts = model_parts

    def GetPythonList(self):
        data = []
        data_mesh = self.solver.model[self.geometry_name]
        data_variable = cs_data_structure.KratosGlobals.GetVariable(self.name)
        for node in data_mesh.Nodes:
            data_value = node.GetSolutionStepValue(data_variable, 0)
            for value in data_value:
                data.append(value)
        return data

    def GetNumpyArray(self):
        # To do
        pass

    def ApplyUpdateToData(self, update):
        data_mesh = self.solver.model[self.geometry_name]
        data_variable = cs_data_structure.KratosGlobals.GetVariable(self.name)
        index = 0
        for node in data_mesh.Nodes:
            updated_value = []
            value = node.GetSolutionStepValue(data_variable, 0)
            for value_i in value:
                updated_value.append(update[index])
                index = index + 1
            node.SetSolutionStepValue(data_variable, 0, updated_value)

    def __add__(self, other):
        pass
