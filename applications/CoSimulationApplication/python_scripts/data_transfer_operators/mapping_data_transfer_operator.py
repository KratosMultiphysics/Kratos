from __future__ import print_function, absolute_import, division  # makes these scripts backward compatible with python 2.6 and 2.7

# Importing the base class
from KratosMultiphysics.CoSimulationApplication.base_classes.co_simulation_base_data_transfer_operator import CoSimulationBaseDataTransferOperator

import KratosMultiphysics.CoSimulationApplication.co_simulation_tools as cs_tools

import KratosMultiphysics
import KratosMultiphysics.MappingApplication as KratosMapping

def Create(settings):
    return MappingDataTransferOperator(settings)

class MappingDataTransferOperator(CoSimulationBaseDataTransferOperator):
    def TransferData(self, from_solver_data, to_solver_data, transfer_options):
        pass

    def PrintInfo(self):
        pass

    def Check(self):
        pass

    @classmethod
    def _GetDefaultSettings(cls):
        this_defaults = cs_tools.cs_data_structure.Parameters("""{
            "mapper_settings" : {
                "mapper_type" : "UNSPECIFIED"
            }
        }""")
        this_defaults.AddMissingParameters(super(MappingDataTransferOperator, cls)._GetDefaultSettings())
        return this_defaults


def GetMapperFlags(transfer_options):
    # currently available mapper-flags
    mapper_flags_dict = {
        "add_values"    : KratosMapping.Mapper.ADD_VALUES,
        "swap_sign"     : KratosMapping.Mapper.SWAP_SIGN,
        "use_transpose" : KratosMapping.Mapper.USE_TRANSPOSE
    }

    mapper_flags = KratosMultiphysics.Flags()
    for i in range(transfer_options.size()):
        flag_name = transfer_options[i].GetString()
        mapper_flags |= mapper_flags_dict[flag_name]

    return mapper_flags


# class KratosIo(CoSimulationBaseIO):
#     def __init__(self, model, custom_settings):
#         self.settings = custom_settings
#         ##settings string in json format
#         default_settings = KratosMultiphysics.Parameters("""{
#             "io_type"  : "kratos_io",
#             "settings" : { }
#         }""")
#         self.settings.ValidateAndAssignDefaults(default_settings)

#         self.mappers = {}

#         super(KratosIo, self).__init__(model, self.settings)
#         ### Constructing the IO for this solver

#     ## ImportCouplingInterfaceData :  used to import data from other solvers
#     #
#     #  @param self            The object pointer.
#     #  @param data_object     python dictionary : configuration of the data to be imported.
#     #                                             data will be imported into this dictionary
#     #  @param from_solver     python object : The solver from which data is to be imported.
#     def ImportCouplingInterfaceData(self, data_object, from_solver):
#         if(from_solver):
#             # exchange data from python cosim solver
#             #if(data_object.Has("mapper_settings")):
#             if data_object.mapper_settings is not None:
#                 origin_geo_name = data_object.origin_data.geometry_name
#                 origin_geo      = from_solver.model[origin_geo_name]
#                 origin_var_name = data_object.origin_data.name

#                 dest_geo_name = data_object.geometry_name
#                 dest_geo      = self.model[dest_geo_name]
#                 dest_var_name = data_object.name

#                 self.__ExecuteMapping(origin_geo_name, origin_geo, origin_var_name, dest_geo_name, dest_geo, dest_var_name, data_object)

#             else:
#                 pass ## We can copy one to one instead of mapping here.
#         else:
#             # import data from remote solver
#             pass

#     ## ExportCouplingInterfaceData :  used to export data to other solvers
#     #
#     #  @param self            The object pointer.
#     #  @param data_object     python dictionary : configuration of the mesh to be exported.
#     #                                             also contains the data to export.
#     #  @param to_solver       python object : The solver to which mesh is to be exported.
#     def ExportCouplingInterfaceData(self, data_object, to_solver):
#         if(to_solver):
#             # exchange data from python cosim solver
#             if data_object.mapper_settings is not None:
#                 dest_geo_name = data_object.destination_data.geometry_name
#                 dest_geo      = to_solver.model[dest_geo_name]
#                 dest_var_name = data_object.destination_data.name

#                 origin_geo_name = data_object.geometry_name
#                 origin_geo      = self.model[origin_geo_name]
#                 origin_var_name = data_object.name

#                 self.__ExecuteMapping(origin_geo_name, origin_geo, origin_var_name, dest_geo_name, dest_geo, dest_var_name, data_object)

#             else:
#                 pass ## We can copy one to one instead of mapping here.
#         else:
#             # import data from remote solver
#             pass

#     def __HasMapper(self, mapper_tuple):
#         return mapper_tuple in self.mappers.keys()

#     def __GetMapper(self, mapper_tuple):
#         return self.mappers[mapper_tuple]

#     def __ExecuteMapping(self, origin_geo_name, origin_geo, origin_var_name, dest_geo_name, dest_geo, dest_var_name, current_data):
#         origin_var = KratosMultiphysics.KratosGlobals.GetVariable(origin_var_name)
#         dest_var = KratosMultiphysics.KratosGlobals.GetVariable(dest_var_name)

#         mapper_flags = GetMapperFlags(current_data)

#         if self.__HasMapper((origin_geo_name,dest_geo_name)):
#             mapper = self.__GetMapper((origin_geo_name,dest_geo_name))
#             mapper.Map(origin_var, dest_var, mapper_flags)
#         elif self.__HasMapper((dest_geo_name,origin_geo_name)):
#             mapper = self.__GetMapper((dest_geo_name,origin_geo_name))
#             mapper.InverseMap(dest_var, origin_var, mapper_flags)
#         else: # Create the mapper
#             # TODO properly get the mapper-settings
#             mapper_settings = KratosMultiphysics.Parameters("""{
#                                                                 "mapper_type" : ""
#                                                             }""")
#             mapper_settings["mapper_type"].SetString(current_data.mapper_settings["type"].GetString())
#             mapper = SetupMapper(origin_geo, dest_geo, mapper_settings)
#             self.mappers[(origin_geo_name,dest_geo_name)] = mapper
#             mapper.Map(origin_var, dest_var, mapper_flags)


# def GetMapperFlags(data_object):
#     mapper_flags_dict = {
#         "add_values"    : MappingApp.Mapper.ADD_VALUES,
#         "swap_sign"     : MappingApp.Mapper.SWAP_SIGN,
#         "use_transpose" : MappingApp.Mapper.USE_TRANSPOSE
#     }

#     mapper_flags = KratosMultiphysics.Flags()
#     if data_object.mapper_settings.Has("flags"):
#         flags = data_object.mapper_settings["flags"]
#         for i in range(flags.size()):
#             flag_name = flags[i].GetString()
#             mapper_flags |= mapper_flags_dict[flag_name]

#     return mapper_flags

# def SetupMapper(geo_origin, geo_destination, mapper_settings):
#     mapper = MappingApp.MapperFactory.CreateMapper(geo_origin, geo_destination, mapper_settings)
#     # check how to create the MPI-version => probably check for IsDistriuted in the Sovler or in the ModelParts (will the external-client solver have an MPIComm?)
#     return mapper
