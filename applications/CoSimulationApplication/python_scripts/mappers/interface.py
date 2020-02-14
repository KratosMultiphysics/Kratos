from KratosMultiphysics.CoSimulationApplication.co_simulation_component import CoSimulationComponent
import KratosMultiphysics.CoSimulationApplication.co_simulation_tools as cs_tools
cs_data_structure = cs_tools.cs_data_structure
import KratosMultiphysics as KM


def Create(parameters):
    return MapperInterface(parameters)


# Class MapperInterface: Interface interpolation with same interpolator type for all modelparts and variables.
class MapperInterface(CoSimulationComponent):
    def __init__(self, parameters):
        super().__init__()

        self.parameters = parameters
        self.settings = parameters["settings"]
        self.mappers = []

    def Initialize(self, interface_from, interface_to):
        super().Initialize()

        # Loop over ModelParts and create mappers
        self.mappers = []
        self.keys = []
        for item_from, item_to in zip(interface_from.model_parts_variables,
                                      interface_to.model_parts_variables):
            key_from = item_from[0]
            key_to = item_to[0]
            self.keys.append((key_from, key_to))  # for PrintInfo

            self.mappers.append(cs_tools.CreateInstance(self.settings))
            self.mappers[-1].Initialize(interface_from.model[key_from],
                                        interface_to.model[key_to])

    def Finalize(self):
        super().Finalize()

        for mapper in self.mappers:
            mapper.Finalize()

    def __call__(self, interface_from, interface_to):
        # Loop over ModelParts and Variables to interpolate
        for i, mapper in enumerate(self.mappers):
            key_from, variables_from = interface_from.model_parts_variables[i]
            key_to, variables_to = interface_to.model_parts_variables[i]
            model_part_from = interface_from.model[key_from]
            model_part_to = interface_to.model[key_to]
            for var_from, var_to in zip(variables_from.list(), variables_to.list()):
                mapper((model_part_from, vars(KM)[var_from.GetString()]),
                       (model_part_to, vars(KM)[var_to.GetString()]))

    def OutputSolutionStep(self):
        for mapper in self.mappers:
            mapper.OutputSolutionStep()

    def PrintInfo(self, indent):
        cs_tools.Print('\t' * indent, "The component ", self.__class__.__name__, " has the following mapper(s):")
        for i, mapper in enumerate(self.mappers):
            mapper.PrintInfo(indent + 1)

            cs_tools.Print('\t' * (indent + 2),
                           f"(which maps ModelPart '{self.keys[i][0]}' to ModelPart '{self.keys[i][1]}')")

