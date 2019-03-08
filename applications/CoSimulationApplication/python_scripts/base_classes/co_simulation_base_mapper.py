# co simulation imports
import co_simulation_tools as tools
# Other imports
import co_simulation_data_structure
cs_data_structure = co_simulation_data_structure.__DATA_STRUCTURE__

##
#  IMPORTANT : This is a BASE CLASS
#               Please do not change anything in this class.
#
#  This class is intended to server as the base class for all the coupled solvers.
class CoSimulationBaseMapper(object):
    def __init__(self, source_geo, destination_geo, settings):
        pass

    ## Map :  Maps the origin_geo to destination_geo
    #
    #  @param self                      The object pointer.
    #  @param from_data                 The origin data to map from
    #  @param to_data                   The destination data to map to
    #  @param mapper_flags              Additional flags for the mapper
    def Map(self, from_data, to_data, mapper_flags):
        pass

    ## InverseMap :  Maps data from destination_geo to the origin_geo
    #
    #  @param self                      The object pointer.
    #  @param from_data                 The origin data to map from
    #  @param to_data                   The destination data to map to
    #  @param mapper_flags              Additional flags for the mapper
    def InverseMap(self, from_data, to_data, mapper_flags):
        pass

    ## UpdateOriginGeometry :  Updates the Origin geometry and re computes the mapping.
    #
    #  @param self                      The object pointer.
    #  @param new_origin_geo            The new geometry with will act as origin geometry.
    def UpdateOriginGeometry(self, new_origin_geo):
        pass

    ## UpdateOriginGeometry :  Updates the Origin geometry and re computes the mapping.
    #
    #  @param self                      The object pointer.
    #  @param new_destination_geo       The new geometry with will act as destination geometry.
    def UpdateDestinationGeometry(self, new_destination_geo):
        pass