import numpy as np

## Class CouplingInterfaceData: Class to hold different properties of the data field contributed in
#                           CoSimulation.
#
class CouplingInterfaceData(object):
    def __init__(self, custom_config, solver):

        default_config = cs_data_structure.Parameters("""
        {
            "name" : "default",
            "dimension" : 0,
            "geometry_name" : "",
            "location_on_mesh":"on_nodes"
        }
        """)
        custom_config.ValidateAndAssignDefaults(default_config)

        self.name = custom_config["name"].GetString()
        self.variable = None
        # self.filters = []
        self.solver = solver
        self.dimension = custom_config["dimension"].GetInt()
        self.location_on_mesh = custom_config["location_on_mesh"].GetString()
        self.mesh_name = custom_config["geometry_name"].GetString()
        self.destination_data = None
        self.origin_data = None
        self.mapper_settings = None

    # def ApplyFilters(self):
    #     for filter in self.filters:
    #         filter.Apply()

    def GetPythonList(self):
        data = []
        data_mesh = self.solver.model[self.mesh_name]
        data_variable = cs_data_structure.KratosGlobals.GetVariable(self.name)
        for node in data_mesh.Nodes:
            data_value = node.GetSolutionStepValue(data_variable,0) #TODO what if non-historical?
            for value in data_value:
                data.append(value)
        return data

    def GetNumpyArray(self):
        return np.array(self.GetPythonList())

    def ApplyUpdateToData(self, update):
        data_mesh = self.solver.model[self.mesh_name]
        data_variable = cs_data_structure.KratosGlobals.GetVariable(self.name)
        index = 0
        for node in data_mesh.Nodes: # #TODO: local nodes to also work in MPI?
            updated_value = []
            value = node.GetSolutionStepValue(data_variable,0)
            # TODO: aditya the data might also be non-historical => GetValue
            for value_i in value:
                updated_value.append(update[index])
                index = index + 1
            node.SetSolutionStepValue(data_variable, 0, updated_value)


