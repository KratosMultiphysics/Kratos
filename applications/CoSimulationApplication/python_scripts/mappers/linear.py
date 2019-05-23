import KratosMultiphysics.CoSimulationApplication.co_simulation_tools as cs_tools
cs_data_structure = cs_tools.cs_data_structure


def Create(parameters):
    return MapperLinear(parameters)


# Class MapperLinear: Linear interpolation.
class MapperLinear(object):
    def __init__(self, _unused):
        super().__init__()
