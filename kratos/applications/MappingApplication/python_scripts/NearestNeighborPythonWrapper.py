# makes KratosMultiphysics backward compatible with python 2.6 and 2.7
from __future__ import print_function, absolute_import, division

from KratosMultiphysics import *
from KratosMultiphysics.MappingApplication import *
# The following check is needed in case of an mpi-parallel compilation that is run serial
try:
    from KratosMultiphysics.mpi import *
except:
    pass

class NearestNeighborPythonWrapper:
    def __init__(self, interface_model_part_origin, interface_model_part_destination):
        self.mapper = NearestNeighborMapper(interface_model_part_origin, interface_model_part_destination)
