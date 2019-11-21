import numpy as np
import math as m
import os.path as path

import KratosMultiphysics as KM
from KratosMultiphysics.CoSimulationApplication.co_simulation_component import CoSimulationComponent
from KratosMultiphysics.CoSimulationApplication.co_simulation_interface import CoSimulationInterface
import KratosMultiphysics.CoSimulationApplication.co_simulation_tools as cs_tools
cs_data_structure = cs_tools.cs_data_structure


def Create(parameters):
    return SolverWrapperAbaqus614(parameters)

class SolverWrapperAbaqus614(CoSimulationComponent):
    def __init__(self, parameters):
        super().__init__()

        #TODO: Continue working on the __init__ function