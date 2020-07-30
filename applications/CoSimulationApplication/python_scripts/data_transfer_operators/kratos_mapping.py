from __future__ import print_function, absolute_import, division  # makes these scripts backward compatible with python 2.6 and 2.7

# Importing the base class
from KratosMultiphysics.CoSimulationApplication.base_classes.co_simulation_data_transfer_operator import CoSimulationDataTransferOperator

# Importing the Kratos Library
import KratosMultiphysics as KM
import KratosMultiphysics.MappingApplication as KratosMapping

# CoSimulation imports
import KratosMultiphysics.CoSimulationApplication.co_simulation_tools as cs_tools
import KratosMultiphysics.CoSimulationApplication.colors as colors

def Create(settings):
    return KratosMappingDataTransferOperator(settings)

class KratosMappingDataTransferOperator(CoSimulationDataTransferOperator):
    """DataTransferOperator that maps values from one interface (ModelPart) to another.
    The mappers of the Kratos-MappingApplication are used
    """
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

    def _ExecuteTransferData(self, from_solver_data, to_solver_data, transfer_options):
        model_part_origin      = from_solver_data.GetModelPart()
        model_part_origin_name = from_solver_data.model_part_name
        variable_origin        = from_solver_data.variable
        identifier_origin      = from_solver_data.solver_name + "." + model_part_origin_name

        model_part_destination      = to_solver_data.GetModelPart()
        model_part_destination_name = to_solver_data.model_part_name
        variable_destination        = to_solver_data.variable
        identifier_destination      = to_solver_data.solver_name + "." + model_part_destination_name

        mapper_flags = self.__GetMapperFlags(transfer_options)
        # TODO in the future automatically add the flags if the values are non-historical

        identifier_tuple         = (identifier_origin, identifier_destination)
        inverse_identifier_tuple = (identifier_destination, identifier_origin)

        if identifier_tuple in self.__mappers:
            self.__mappers[identifier_tuple].Map(variable_origin, variable_destination, mapper_flags)
        elif inverse_identifier_tuple in self.__mappers:
            self.__mappers[inverse_identifier_tuple].InverseMap(variable_destination, variable_origin, mapper_flags)
        else:
            if model_part_origin.IsDistributed() or model_part_destination.IsDistributed():
                mapper_create_fct = KratosMapping.MapperFactory.CreateMPIMapper
            else:
                mapper_create_fct = KratosMapping.MapperFactory.CreateMapper

            if self.echo_level > 0:
                info_msg  = "Creating Mapper:\n"
                info_msg += '    Origin: ModelPart "{}" of solver "{}"\n'.format(model_part_origin_name, from_solver_data.solver_name)
                info_msg += '    Destination: ModelPart "{}" of solver "{}"'.format(model_part_destination_name, to_solver_data.solver_name)

                cs_tools.cs_print_info(colors.bold(self._ClassName()), info_msg)

            self.__mappers[identifier_tuple] = mapper_create_fct(model_part_origin, model_part_destination, self.settings["mapper_settings"].Clone()) # Clone is necessary because the settings are validated and defaults assigned, which could influence the creation of other mappers
            self.__mappers[identifier_tuple].Map(variable_origin, variable_destination, mapper_flags)

    def _Check(self, from_solver_data, to_solver_data):
        def CheckData(data_to_check):
            if data_to_check.location != "node_historical":
                raise Exception('Currently only historical nodal values are supported by the "{}"\nChecking ModelPart "{}" of solver "{}"'.format(self._ClassName(), data_to_check.model_part_name, data_to_check.solver_name))

        CheckData(from_solver_data)
        CheckData(to_solver_data)

        # TODO in the future also non-historical nodal values will be supported, but this still requires some improvements in the MappingApp

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
