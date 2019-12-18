import numpy as np

import KratosMultiphysics.CoSimulationApplication.co_simulation_tools as cs_tools
from KratosMultiphysics.CoSimulationApplication.mappers.interpolator import MapperInterpolator
cs_data_structure = cs_tools.cs_data_structure

import time
from contextlib import contextmanager
@contextmanager
def timer(name=None, t=0, n=0, ms=False):
    startTime = time.time()
    yield
    elapsedTime = time.time() - startTime
    if ms:
        s = '\n' * n + '\t' * t + f'{elapsedTime * 1000:.2f}ms'
        s.replace(',', ' ')
    else:
        s = '\n' * n + '\t' * t + f'{elapsedTime:.1f}s'
    if name is not None:
        s += f' - {name}'
    s += '\n' * n
    print(s)

def Create(parameters):
    return MapperNearest(parameters)


# Class MapperNearest: 3D nearest-neighbour interpolation.
class MapperNearest(MapperInterpolator):
    def __init__(self, parameters):
        super().__init__(parameters)

        self.n_nearest = 1

    def Initialize(self, model_part_from, model_part_to):
        super().Initialize(model_part_from, model_part_to)

        self.nearest = self.nearest.reshape(-1, 1)
        self.coeffs = np.ones((self.n_to, 1))
