# co simulation imports
import co_simulation_tools as tools
# Importing the CoSimulation application
from CoSimulationApplication import *
# Other imports
import co_simulation_data_structure
cs_data_structure = co_simulation_data_structure.__KRATOS_DATA_STRUCTURE__
import collections

##
#  IMPORTANT : This is a BASE CLASS
#               Please do not change anything in this class.
#
#  This class is intended to server as the base class for all the coupled solvers.
class CoSimulationBaseMapper(object):
    def __init__(self, source_geo, destination_geo, settings):
        pass

    ## Map :  Maps the from_data to to_data
    #
    #  @param self                      The object pointer.
    #  @param from_data                 The origin data to map from
    #  @param to_data                   The destination data to map to
    #  @param mapper_flags              Additional flags for the mapper
    def Map(self, from_data, to_data, mapper_flags):
        pass

    def UpdateSourceGeometry(self):
        pass

    def UpdateDestinationGeometry(self):
        pass
