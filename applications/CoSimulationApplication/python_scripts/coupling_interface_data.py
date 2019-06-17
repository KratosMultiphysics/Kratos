
# Other imports
from . import co_simulation_tools as cs_tools
import numpy as np

## Class CouplingInterfaceData: Class to hold different properties of the data field contributed in
#                           CoSimulation.
#
class CouplingInterfaceData(object):
    def __init__(self, custom_config, model):

        default_config = cs_tools.cs_data_structure.Parameters("""{
            "model_part_name" : "UNSPECIFIED",
            "variable_name"   : "UNSPECIFIED",
            "location"        : "node_historical",
            "dimension"       : -1
        }""")
        custom_config.ValidateAndAssignDefaults(default_config)

        self.variable = cs_tools.cs_data_structure.KratosGlobals.GetVariable(custom_config["variable_name"].GetString())
        self.model = model
        self.dimension = custom_config["dimension"].GetInt() # TODO check that sth was assigned
        self.location = custom_config["location"].GetString()
        self.model_part_name = custom_config["model_part_name"].GetString()
        # TODO check DOMAIN_SIZE against the dimension

    def GetModelPart(self):
        return self.model[self.model_part_name]

    def GetPythonList(self, solution_step_index=0):
        model_part = self.GetModelPart()
        data = [0]*len(model_part.Nodes)*self.dimension
        node_index = 0
        for node in model_part.Nodes:
            data_value = node.GetSolutionStepValue(self.variable,solution_step_index) #TODO what if non-historical?
            for i in range(self.dimension):
                data[node_index*self.dimension + i] = data_value[i]
            node_index+=1
        return data

    def GetNumpyArray(self, solution_step_index=0):
        return np.asarray(self.GetPythonList(solution_step_index), dtype=np.float64)

    def ApplyUpdateToData(self, update):
        model_part = self.GetModelPart()
        node_index = 0
        updated_value = [0]*3 # TODO this is a hack, find better solution! => check var_type
        for node in model_part.Nodes: # #TODO: local nodes to also work in MPI?
            # TODO: aditya the data might also be non-historical => GetValue
            for i in range(self.dimension):
                updated_value[i] = update[node_index*self.dimension + i]

            node.SetSolutionStepValue(self.variable, 0, updated_value)
            node_index += 1
