# Importing the base class
from KratosMultiphysics.CoSimulationApplication.base_classes.co_simulation_data_transfer_operator import CoSimulationDataTransferOperator

# Importing the Kratos Library
import KratosMultiphysics as KM
import KratosMultiphysics.MappingApplication as KratosMapping

# CoSimulation imports
import KratosMultiphysics.CoSimulationApplication.co_simulation_tools as cs_tools
import KratosMultiphysics.CoSimulationApplication.colors as colors

# other imports
from KratosMultiphysics.CoSimulationApplication.utilities import data_communicator_utilities
from time import time

def Create(*args):
    return Kratos3D1DDataTransferOperator(*args)

class Kratos3D1DDataTransferOperator(CoSimulationDataTransferOperator):
    """DataTransferOperator that maps values from one interface (ModelPart) to another.
    The mappers of the Kratos-MappingApplication are used
    """

    # initializing the static members necessary for MPI
    # initializing on the fly does not work and leads to memory problems
    # as the members are not proberly saved and randomly destucted!
    __dummy_model = None
    __rank_zero_model_part = None

    def __init__(self, settings, parent_coupled_solver_data_communicator):
        if not settings.Has("mapper_settings"):
            raise Exception('No "mapper_settings" provided!')
        super().__init__(settings, parent_coupled_solver_data_communicator)
        self.__mappers = {}

    def _ExecuteTransferData(self, from_solver_data, to_solver_data, transfer_options):
        model_part_origin_name = from_solver_data.model_part_name
        variable_origin        = from_solver_data.variable
        identifier_origin      = from_solver_data.solver_name + "." + model_part_origin_name

        model_part_destination_name = to_solver_data.model_part_name
        variable_destination        = to_solver_data.variable
        identifier_destination      = to_solver_data.solver_name + "." + model_part_destination_name

        mapper_flags = self.__GetMapperFlags(transfer_options, from_solver_data, to_solver_data)

        identifier_tuple         = (identifier_origin, identifier_destination)
        inverse_identifier_tuple = (identifier_destination, identifier_origin)

        if identifier_tuple in self.__mappers:
            self.__mappers[identifier_tuple].Map(variable_origin, variable_destination, mapper_flags)
        elif inverse_identifier_tuple in self.__mappers:
            self.__mappers[inverse_identifier_tuple].InverseMap(variable_destination, variable_origin, mapper_flags)
        else:
            model_part_origin      = self.__GetModelPartFromInterfaceData(from_solver_data)
            model_part_destination = self.__GetModelPartFromInterfaceData(to_solver_data)

            if model_part_origin.IsDistributed() or model_part_destination.IsDistributed():
                mapper_create_fct = python_mapper_factory.CreateMPIMapper
            else:
                mapper_create_fct = python_mapper_factory.CreateMapper

            if self.echo_level > 0:
                info_msg  = "Creating Mapper:\n"
                info_msg += '    Origin: ModelPart "{}" of solver "{}"\n'.format(model_part_origin_name, from_solver_data.solver_name)
                info_msg += '    Destination: ModelPart "{}" of solver "{}"'.format(model_part_destination_name, to_solver_data.solver_name)

                cs_tools.cs_print_info(colors.bold(self._ClassName()), info_msg)

            mapper_creation_start_time = time()
            self.__mappers[identifier_tuple] = mapper_create_fct(model_part_origin, model_part_destination, self.settings["mapper_settings"].Clone()) # Clone is necessary because the settings are validated and defaults assigned, which could influence the creation of other mappers

            if self.echo_level > 2:
                cs_tools.cs_print_info(colors.bold(self._ClassName()), "Creating Mapper took: {0:.{1}f} [s]".format(time()-mapper_creation_start_time,2))
            self.__mappers[identifier_tuple].Map(variable_origin, variable_destination, mapper_flags)

    def _Check(self, from_solver_data, to_solver_data):
        def CheckData(data_to_check):
            if "node" not in data_to_check.location:
                raise Exception('Mapping only supports nodal values!"{}"\nChecking ModelPart "{}" of solver "{}"'.format(self._ClassName(), data_to_check.model_part_name, data_to_check.solver_name))

        CheckData(from_solver_data)
        CheckData(to_solver_data)

    @classmethod
    def _GetDefaultParameters(cls):
        this_defaults = KM.Parameters("""{
            "mapper_settings" : {
                "mapper_type" : "UNSPECIFIED"
            }
        }""")
        this_defaults.AddMissingParameters(super()._GetDefaultParameters())
        return this_defaults

    @classmethod
    def _GetListAvailableTransferOptions(cls):
        return ["swap_sign"]

    @staticmethod
    def __GetModelPartFromInterfaceData(interface_data):
        """If the solver does not exist on this rank, then pass a
        dummy ModelPart to the Mapper that has a DataCommunicator
        that is not defined on this rank
        """
        if interface_data.IsDefinedOnThisRank():
            return interface_data.GetModelPart()
        else:
            if KM.IsDistributedRun():
                raise Exception("This function can only be called when Kratos is running in MPI!")
