import KratosMultiphysics.CoSimulationApplication.co_simulation_tools as cs_tools
cs_data_structure = cs_tools.cs_data_structure


def Create(parameters):
    return MapperInterface(parameters)


# Class MapperInterface: Interface interpolation with same interpolator type for all modelparts.
class MapperInterface(object):
    def __init__(self, parameters):
        super().__init__()

        self.parameters = parameters
        self.settings = parameters["settings"]

	def Initialize(self, interface_from, interface_to):
		# Loop over modelparts and create mappers
		keys_from = [_[0] for _ in interface_from.model_parts_variables]
		keys_to = [_[0] for _ in interface_to.model_parts_variables]
		for i in range(len(keys_from)):
			self.mappers[i] = cs_tools.CreateInstance(self.settings)
			self.mappers[i].Initialize(interface_from.model[keys_from[i]], interface_to.model[keys_to[i]])
			
	def Finalize(self):
		for mapper in self.mappers:
			mapper.Finalize()

	def __call__(self, interface_from, interface_to):
		# Loop over modelparts and variables to interpolate
		keys_from = [_[0] for _ in interface_from.model_parts_variables]
		keys_to = [_[0] for _ in interface_to.model_parts_variables]
		variables_from = [_[1].list() for _ in interface_from.model_parts_variables]
		variables_to = [_[1].list() for _ in interface_to.model_parts_variables]
		
