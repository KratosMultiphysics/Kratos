from __future__ import print_function, absolute_import, division  # makes these scripts backward compatible with python 2.6 and 2.7

# Importing the base class
from KratosMultiphysics.CoSimulationApplication.base_classes.co_simulation_data_transfer_operator import CoSimulationDataTransferOperator

# Importing the Kratos Library
import KratosMultiphysics as KM
import KratosMultiphysics.MappingApplication as KratosMapping

def Create(settings):
    return KratosMappingDataTransferOperator(settings)

class KratosMappingDataTransferOperator(CoSimulationDataTransferOperator):

    # currently available mapper-flags aka transfer-options
    __mapper_flags_dict = {
        "add_values"    : KratosMapping.Mapper.ADD_VALUES,
        "swap_sign"     : KratosMapping.Mapper.SWAP_SIGN,
        "use_transpose" : KratosMapping.Mapper.USE_TRANSPOSE
    }

    def __init__(self, settings):
        if not settings.Has("mapper_settings"):
            raise Exception('No "mapper_settings" provided!')
        super(KratosMappingDataTransferOperator, self).__init__(settings)
        self.__mappers = {}

    def TransferData(self, from_solver_data, to_solver_data, transfer_options):
        # TODO check location of data => should coincide with the one for the mapper
        # or throw if it is not in a suitable location (e.g. on the ProcessInfo)

        self._CheckAvailabilityTransferOptions(transfer_options)

        model_part_origin      = from_solver_data.GetModelPart()
        model_part_origin_name = model_part_origin.Name
        variable_origin        = from_solver_data.variable

        model_part_destinatinon      = to_solver_data.GetModelPart()
        model_part_destinatinon_name = model_part_destinatinon.Name
        variable_destination         = to_solver_data.variable

        mapper_flags = self.__GetMapperFlags(transfer_options)

        name_tuple         = (model_part_origin_name, model_part_destinatinon_name)
        inverse_name_tuple = (model_part_destinatinon_name, model_part_origin_name)

        if name_tuple in self.__mappers:
            self.__mappers[name_tuple].Map(variable_origin, variable_destination, mapper_flags)
        elif inverse_name_tuple in self.__mappers:
            self.__mappers[inverse_name_tuple].InverseMap(variable_destination, variable_origin, mapper_flags)
        else:
            if model_part_origin.IsDistributed() or model_part_destinatinon.IsDistributed():
                mapper_create_fct = KratosMapping.MapperFactory.CreateMPIMapper
            else:
                mapper_create_fct = KratosMapping.MapperFactory.CreateMapper

            self.__mappers[name_tuple] = mapper_create_fct(model_part_origin, model_part_destinatinon, self.settings["mapper_settings"].Clone()) # Clone is necessary here bcs settings are influenced among mappers otherwise. TODO check in the MapperFactory how to solve this better
            self.__mappers[name_tuple].Map(variable_origin, variable_destination, mapper_flags)

    @classmethod
    def _GetDefaultSettings(cls):
        this_defaults = KM.Parameters("""{
            "mapper_settings" : {
                "mapper_type" : "UNSPECIFIED"
            }
        }""")
        this_defaults.AddMissingParameters(super(KratosMappingDataTransferOperator, cls)._GetDefaultSettings())
        return this_defaults

    @classmethod
    def _GetListAvailableTransferOptions(cls):
        return cls.__mapper_flags_dict.keys()

    def __GetMapperFlags(self, transfer_options):
        mapper_flags = KM.Flags()
        for flag_name in transfer_options.GetStringArray():
            mapper_flags |= self.__mapper_flags_dict[flag_name]

        return mapper_flags
