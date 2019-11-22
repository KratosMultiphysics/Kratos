import KratosMultiphysics.CoSimulationApplication.co_simulation_tools as cs_tools
cs_data_structure = cs_tools.cs_data_structure


def Create(parameters):
    return MapperLinear(parameters)


# Class MapperLinear: Linear interpolation.
class MapperLinear(object):
    def __init__(self, _unused):
        super().__init__()

    def Initialize(self, model_part_from, model_part_to):
        raise NotImplementedError

    def Finalize(self):
        pass

    def __call__(self, args_from, args_to):
        model_part_from, var_from = args_from
        model_part_to, var_to = args_to
        raise NotImplementedError
