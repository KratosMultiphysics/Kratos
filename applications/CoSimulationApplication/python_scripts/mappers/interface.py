import KratosMultiphysics.CoSimulationApplication.co_simulation_tools as cs_tools
cs_data_structure = cs_tools.cs_data_structure


def Create(parameters):
    return MapperInterface(parameters)


# Class MapperInterface: Interface interpolation with same interpolator type for all modelparts.
class MapperInterface(object):
    def __init__(self, _unused):
        super().__init__()
