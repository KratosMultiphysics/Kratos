from KratosMultiphysics import *
from KratosMultiphysics.SwimmingDEMApplication import *
from KratosMultiphysics.DEMApplication import *
import ethier_benchmark_algorithm
import swimming_DEM_procedures as SDP
import math
import numpy as np
BaseAlgorithm = ethier_benchmark_algorithm.Algorithm

class Algorithm(BaseAlgorithm):
    def __init__(self, varying_parameters = Parameters("{}")):
        BaseAlgorithm.__init__(self, varying_parameters)

    def GetFieldUtility(self):
        self.flow_field = ProductOfSines()
        space_time_set = SpaceTimeSet()
        self.field_utility = FluidFieldUtility(space_time_set, self.flow_field, 1000.0, 1e-6)
        return self.field_utility