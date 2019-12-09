from scipy.spatial import cKDTree

import numpy as np

import KratosMultiphysics.CoSimulationApplication.co_simulation_tools as cs_tools
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
    return MapperCombined(parameters)


# Class MapperCombined
class MapperCombined(object):
    def __init__(self, parameters):
        """
        TODO - add some explanation.
        """
        super().__init__()

        self.settings = parameters['settings']

        # create all mappers
        self.mappers = []
        for par in self.settings['mappers'].list():
            self.mappers.append(cs_tools.CreateInstance(par))

        # check that exactly one mapper is an interpolator
        counter = 0
        for i, mapper in enumerate(self.mappers):
            if mapper.interpolator:
                self.index = i
                counter += 1
        if counter != 1:
            raise ValueError(f'{counter} interpolators found instead of 1')
        # *** add this attribute to every mapper!

    def Initialize(self, model_part_from, model_part_to):
        # initialize usptream transformers
        mps_from = []

        # *** TODO

        # initialize downstream transformers
        mps_to = []

        # *** TODO

        mps_to.reverse()
        self.model_parts = mps_from + mps_to

        # initialize interpolator
        self.mappers[self.index].Initialize(mps_from[-1], mps_to[0])

    def Finalize(self):
        pass

    def __call__(self, args_from, args_to):
        mps_from = [args_from[0]] + self.model_parts
        mps_to = self.model_parts + [args_to[0]]

        for i, mapper in enumerate(self.mappers):
            var_from = args_from[1]
            var_to = args_to[1]
            if i > self.index:
                var_from = args_to[1]
            if i < self.index:
                var_to = args_from[1]
            mapper((mps_from[i], var_from), (mps_to[i], var_to))
