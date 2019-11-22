import KratosMultiphysics.CoSimulationApplication.co_simulation_tools as cs_tools
cs_data_structure = cs_tools.cs_data_structure
import KratosMultiphysics as KM


def Create(parameters):
    return MapperInterface(parameters)


# Class MapperInterface: Interface interpolation with same interpolator type for all modelparts and variables.
class MapperInterface(object):
    def __init__(self, parameters):
        super().__init__()

        self.parameters = parameters
        self.settings = parameters["settings"]

    def Initialize(self, interface_from, interface_to):
        # Loop over ModelParts and create mappers
        self.mappers = []
        for item_from, item_to in zip(interface_from.model_parts_variables,
                                      interface_to.model_parts_variables):
            key_from, = item_from
            key_to, = item_to
            self.mappers.append(cs_tools.CreateInstance(self.settings))
            self.mappers[-1].Initialize(interface_from.model[key_from],
                                        interface_to.model[key_to])

    def Finalize(self):
        for mapper in self.mappers:
            mapper.Finalize()

    def __call__(self, interface_from, interface_to):
        # Loop over modelparts and variables to interpolate
        for i, mapper in enumerate(self.mappers):
            key_from, variables_from = interface_from.model_parts_variables[i]
            key_to, variables_to = interface_to.model_parts_variables[i]
            model_part_from = interface_from.model[key_from]
            model_part_to = interface_from.model[key_to]
            for var_from, var_to in zip(variables_from, variables_to):
                mapper((model_part_from, vars(KM)[var_from.GetString()]),
                       (model_part_to, vars(KM)[var_to.GetString()]))
