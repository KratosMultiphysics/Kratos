import numpy as np

import KratosMultiphysics.CoSimulationApplication.co_simulation_tools as cs_tools
from KratosMultiphysics.CoSimulationApplication.mappers.interpolator import MapperInterpolator
cs_data_structure = cs_tools.cs_data_structure


def Create(parameters):
    return MapperNearest(parameters)


# Class MapperNearest: nearest-neighbour interpolation.
class MapperNearest(MapperInterpolator):
    def __init__(self, parameters):
        super().__init__(parameters)

        self.n_nearest = 1

    def Initialize(self, model_part_from, model_part_to):
        super().Initialize(model_part_from, model_part_to)

        # self.nearest = self.nearest.reshape(-1, 1)
        self.coeffs = np.ones((self.n_to, 1))
