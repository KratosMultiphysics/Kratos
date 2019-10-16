from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7

# Python factory for all the currently available Mappers in Kratos
# The intention is to give the users a unique place to create Mappers
# The goal is to implement the Mappers from the other Apps also in the
# MappingApp (which inherently also work in MPI) and replace them in
# the long run.
# This way users won't notice / won't have to change their code

import KratosMultiphysics as KM
from KratosMultiphysics.MappingApplication import MapperFactory

from importlib import import_module

available_mappers = {
    "empire_nearest_neighbor" : "empire_mapper_wrapper",
    "empire_nearest_element"  : "empire_mapper_wrapper",
    "empire_barycentric"      : "empire_mapper_wrapper",
    "empire_dual_mortar"      : "empire_mapper_wrapper",
    "empire_mortar"           : "empire_mapper_wrapper"
}

def CreateMapper(model_part_origin, model_part_destination, mapper_settings):
    if not isinstance(mapper_settings, KM.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")

    mapper_type = mapper_settings["mapper_type"].GetString()

    is_distributed = model_part_origin.IsDistributed() or model_part_destination.IsDistributed()

    if is_distributed:
        list_avail_mappers = MapperFactory.GetRegisteredMPIMapperNames()
        dist_info = "MPI"
    else:
        list_avail_mappers = MapperFactory.GetRegisteredMapperNames()
        dist_info = "serial (non-MPI)"

    # use the MappingApp if it has the requested mapper
    if is_distributed:
        if MapperFactory.HasMPIMapper(mapper_type):
            return MapperFactory.CreateMPIMapper(model_part_origin, model_part_destination, mapper_settings)
    elif MapperFactory.HasMapper(mapper_type):
        return MapperFactory.CreateMapper(model_part_origin, model_part_destination, mapper_settings)

    if mapper_type in available_mappers:
        mapper_module_name = available_mappers[mapper_type]

        if mapper_type.startswith("empire"):
            raw_mapper_name = mapper_type[7:]
            if raw_mapper_name in list_avail_mappers:
                info_msg  = "Using a mapper from Empire ({}), even though the same mapper is available in Kratos\n".format(raw_mapper_name)
                info_msg += "Consider switching for improved performance, less overead and native MPI support"
                KM.Logger.PrintInfo("Python-Mapper-Factory", info_msg)

        return import_module("KratosMultiphysics.MappingApplication." + mapper_module_name).CreateMapper(model_part_origin, model_part_destination, mapper_settings)
    else:
        list_avail_mappers.extend(available_mappers.keys())

        err_msg  = 'The requested mapper "{}" is not available in {}\n'.format(mapper_type, dist_info)
        err_msg += 'The following mappers are available:'
        for avail_mapper in list_avail_mappers:
            err_msg += '\n\t{}'.format(avail_mapper)
        raise Exception(err_msg)
