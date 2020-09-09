# Python factory for all the currently available Mappers in Kratos
# The intention is to give the users a unique place to create Mappers
# The goal is to implement the Mappers from the other Apps also in the
# MappingApp (which inherently also work in MPI) and replace them in
# the long run.
# This way users won't notice / won't have to change their code

import KratosMultiphysics as KM
from KratosMultiphysics.MappingApplication import MapperFactory

def CreateMapper(model_part_origin, model_part_destination, mapper_settings):
    if not isinstance(mapper_settings, KM.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")

    mapper_type = mapper_settings["mapper_type"].GetString()

    # use the MappingApp if it has the requested mapper
    if MapperFactory.HasMapper(mapper_type):
        return MapperFactory.CreateMapper(model_part_origin, model_part_destination, mapper_settings)

    list_avail_mappers = MapperFactory.GetRegisteredMapperNames()

    err_msg  = 'The requested mapper "{}" is not available in serial (non-MPI)\n'.format(mapper_type)
    err_msg += 'The following mappers are available:'
    for avail_mapper in list_avail_mappers:
        err_msg += '\n\t{}'.format(avail_mapper)
    raise Exception(err_msg)


def CreateMPIMapper(model_part_origin, model_part_destination, mapper_settings):
    if not isinstance(mapper_settings, KM.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")

    mapper_type = mapper_settings["mapper_type"].GetString()

    # use the MappingApp if it has the requested mapper
    if MapperFactory.HasMPIMapper(mapper_type):
        return MapperFactory.CreateMPIMapper(model_part_origin, model_part_destination, mapper_settings)

    list_avail_mappers = MapperFactory.GetRegisteredMPIMapperNames()

    err_msg  = 'The requested mapper "{}" is not available in MPI\n'.format(mapper_type)
    err_msg += 'The following mappers are available:'
    for avail_mapper in list_avail_mappers:
        err_msg += '\n\t{}'.format(avail_mapper)
    raise Exception(err_msg)
