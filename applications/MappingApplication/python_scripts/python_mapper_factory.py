from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7

# Python factory for all the currently available Mappers in Kratos
# The intention is to give the users a unique place to create Mappers
# The goal is to implement the Mappers from the other Apps also in the
# MappingApp (which inherently also work in MPI) and replace them in
# the long run.
# This way users won't notice / won't have to change their code

import KratosMultiphysics
from KratosMultiphysics.MappingApplication import MapperFactory

def _CreateCoreMortarMapper(model_part_origin, model_part_destination, mapper_settings):
    # return core mortar mapper
    raise NotImplementedError

def _CreateApproximateMortarMapper(model_part_origin, model_part_destination, mapper_settings):
    # return mapper from FSIApp
    raise NotImplementedError

def _CreateVertexMorphingMapper(model_part_origin, model_part_destination, mapper_settings):
    # return mapper from ShapeOptApp
    raise NotImplementedError

def _CreateEmpireMortarMapper(model_part_origin, model_part_destination, mapper_settings):
    # return mortar mapper from Empire
    raise NotImplementedError

available_mappers = {
    "mortar"                      : _CreateCoreMortarMapper,
    "approximate_mortar"          : _CreateApproximateMortarMapper,
    "vertex_morphing"             : _CreateVertexMorphingMapper,
    "empire_mortar"               : _CreateEmpireMortarMapper
}

def CreateMapper(model_part_origin, model_part_destination, mapper_settings):
    if not isinstance(mapper_settings, KratosMultiphysics.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")

    mapper_type = mapper_settings["mapper_type"].GetString()

    if MapperFactory.HasMapper(mapper_type): # use the MappingApp if it has the requested mapper
        return MapperFactory.CreateMapper(model_part_origin, model_part_destination, mapper_settings)
    elif mapper_type in available_mappers:
        return available_mappers[mapper_type](model_part_origin, model_part_destination, mapper_settings)
    else:
        list_mappers = MapperFactory.GetRegisteredMapperNames()
        list_mappers.extend(available_mappers.keys())

        err_msg  = 'The requested mapper "{}" is not available\n'.format(mapper_type)
        err_msg += 'The following mappers are available:'
        for avail_mapper in list_mappers:
            err_msg += '\n\t{}'.format(avail_mapper)
        raise Exception(err_msg)

def CreateMPIMapper(model_part_origin, model_part_destination, mapper_settings):
    return MapperFactory.CreateMPIMapper(model_part_origin, model_part_destination, mapper_settings)
