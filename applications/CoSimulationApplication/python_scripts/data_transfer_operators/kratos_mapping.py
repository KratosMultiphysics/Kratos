# Importing the base class
from KratosMultiphysics.CoSimulationApplication.base_classes.co_simulation_data_transfer_operator import CoSimulationDataTransferOperator

# Importing the Kratos Library
import KratosMultiphysics as KM
from KratosMultiphysics.MappingApplication import python_mapper_factory

# CoSimulation imports
import KratosMultiphysics.CoSimulationApplication.co_simulation_tools as cs_tools
import KratosMultiphysics.CoSimulationApplication.colors as colors

# Other imports
from KratosMultiphysics.CoSimulationApplication.utilities import data_communicator_utilities
from time import time
from dataclasses import dataclass
from typing import Union, Tuple

@dataclass
class SolverData:
    model_part_origin: KM.ModelPart
    model_part_destination: KM.ModelPart
    model_part_origin_name: str
    model_part_destination_name: str
    variable_origin: Union[KM.DoubleVariable, KM.Array1DVariable3]
    variable_destination: Union[KM.DoubleVariable, KM.Array1DVariable3]
    mapper_flags: KM.Flags
    identifier_origin: str
    identifier_destination: str
    identifier_tuple: Tuple[str, str]
    inverse_identifier_tuple: Tuple[str, str]

def Create(*args):
    return KratosMappingDataTransferOperator(*args)
class KratosMappingDataTransferOperator(CoSimulationDataTransferOperator):
    """DataTransferOperator that maps values from one interface (ModelPart) to another.
    The mappers of the Kratos-MappingApplication are used
    """
    # currently available mapper-flags aka transfer-options
    __mapper_flags_dict = {
        "add_values"    : KM.Mapper.ADD_VALUES,
        "swap_sign"     : KM.Mapper.SWAP_SIGN,
        "use_transpose" : KM.Mapper.USE_TRANSPOSE
    }

    # initializing the static members necessary for MPI
    # initializing on the fly does not work and leads to memory problems
    # as the members are not proberly saved and randomly destucted!
    __dummy_model = None
    __rank_zero_model_part = None

    def __init__(self, settings, parent_coupled_solver_data_communicator):
        if not settings.Has("mapper_settings"):
            raise Exception('No "mapper_settings" provided!')
        super().__init__(settings, parent_coupled_solver_data_communicator)
        self._mappers = {}

    def _ExecuteTransferData(self, from_solver_data, to_solver_data, transfer_options):
        solver_data = self._PrepareSolverData(from_solver_data, to_solver_data, transfer_options)

        if solver_data.identifier_tuple in self._mappers:
            self._mappers[solver_data.identifier_tuple].Map(solver_data.variable_origin, solver_data.variable_destination, solver_data.mapper_flags)
        elif solver_data.inverse_identifier_tuple in self._mappers:
            self._mappers[solver_data.inverse_identifier_tuple].InverseMap(solver_data.variable_destination, solver_data.variable_origin, solver_data.mapper_flags)
        else:
            model_part_origin      = self._GetModelPartFromInterfaceData(from_solver_data)
            model_part_destination = self._GetModelPartFromInterfaceData(to_solver_data)

            mapper_create_fct = self._DefineMapperFunction(model_part_origin, model_part_destination)

            if self.echo_level > 0:
                info_msg  = "Creating Mapper:\n"
                info_msg += '    Origin: ModelPart "{}" of solver "{}"\n'.format(solver_data.model_part_origin_name, from_solver_data.solver_name)
                info_msg += '    Destination: ModelPart "{}" of solver "{}"'.format(solver_data.model_part_destination_name, to_solver_data.solver_name)

                cs_tools.cs_print_info(colors.bold(self._ClassName()), info_msg)

            mapper_creation_start_time = time()
            self._mappers[solver_data.identifier_tuple] = mapper_create_fct(model_part_origin, model_part_destination, self.settings["mapper_settings"].Clone()) # Clone is necessary because the settings are validated and defaults assigned, which could influence the creation of other mappers

            if self.echo_level > 2:
                cs_tools.cs_print_info(colors.bold(self._ClassName()), "Creating Mapper took: {0:.{1}f} [s]".format(time()-mapper_creation_start_time,2))
            self._mappers[solver_data.identifier_tuple].Map(solver_data.variable_origin, solver_data.variable_destination, solver_data.mapper_flags)

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
        return cls.__mapper_flags_dict.keys()

    def _DefineMapperFunction(self, model_part_origin, model_part_destination):
        """
        Define the mapper function to be used

        Args:
            self: The instance of the class.
            model_part_origin (ModelPart): The model part to transfer from.
            model_part_destination (ModelPart): The model part to transfer to.

        Returns:
            function: The function to be used for the mapping
        """
        if model_part_origin.IsDistributed() or model_part_destination.IsDistributed():
            return python_mapper_factory.CreateMPIMapper
        else:
            return python_mapper_factory.CreateMapper

    def _PrepareSolverData(self, from_solver_data, to_solver_data, transfer_options):
        """
        Prepare the solver data

        Args:
            self: The instance of the class.
            from_solver_data (CoSimulationData): The data from the solver to transfer from.
            to_solver_data (CoSimulationData): The data from the solver to transfer to.
            transfer_options (KM.Flags): The flags to be used for the mapping.

        Returns:
            SolverData: Named tuple containing all the extracted values.
        """
        # Get the from solver data
        model_part_origin_name = from_solver_data.model_part_name
        variable_origin        = from_solver_data.variable
        identifier_origin      = from_solver_data.solver_name + "." + model_part_origin_name

        # Get the to solver data
        model_part_destination_name = to_solver_data.model_part_name
        variable_destination        = to_solver_data.variable
        identifier_destination      = to_solver_data.solver_name + "." + model_part_destination_name

        # Get the mapper flags
        mapper_flags = self.__GetMapperFlags(transfer_options, from_solver_data, to_solver_data)

        # Get the identifier tuple
        identifier_tuple         = (identifier_origin, identifier_destination)
        inverse_identifier_tuple = (identifier_destination, identifier_origin)

        # Get the model parts
        model_part_origin      = self._GetModelPartFromInterfaceData(from_solver_data)
        model_part_destination = self._GetModelPartFromInterfaceData(to_solver_data)

        # Return the solver data as data class
        return SolverData(
            model_part_origin=model_part_origin,
            model_part_destination=model_part_destination,
            model_part_origin_name=model_part_origin_name,
            model_part_destination_name=model_part_destination_name,
            variable_origin=variable_origin,
            variable_destination=variable_destination,
            mapper_flags=mapper_flags,
            identifier_origin=identifier_origin,
            identifier_destination=identifier_destination,
            identifier_tuple=identifier_tuple,
            inverse_identifier_tuple=inverse_identifier_tuple
        )

    @staticmethod
    def _GetModelPartFromInterfaceData(interface_data):
        """If the solver does not exist on this rank, then pass a
        dummy ModelPart to the Mapper that has a DataCommunicator
        that is not defined on this rank
        """
        if interface_data.IsDefinedOnThisRank():
            return interface_data.GetModelPart()
        else:
            return KratosMappingDataTransferOperator.__GetRankZeroModelPart()

    def __GetMapperFlags(self, transfer_options, from_solver_data, to_solver_data):
        mapper_flags = KM.Flags()
        for flag_name in transfer_options.GetStringArray():
            mapper_flags |= self.__mapper_flags_dict[flag_name]
        if from_solver_data.location == "node_non_historical":
            mapper_flags |= KM.Mapper.FROM_NON_HISTORICAL
        if to_solver_data.location == "node_non_historical":
            mapper_flags |= KM.Mapper.TO_NON_HISTORICAL

        return mapper_flags

    @staticmethod
    def __GetRankZeroModelPart():
        if not KM.IsDistributedRun():
            raise Exception("This function can only be called when Kratos is running in MPI!")

        if KratosMappingDataTransferOperator.__rank_zero_model_part is None:
            model = KM.Model()
            rank_zero_model_part = model.CreateModelPart("rank_zero")

            from KratosMultiphysics.mpi import ModelPartCommunicatorUtilities
            ModelPartCommunicatorUtilities.SetMPICommunicator(rank_zero_model_part, data_communicator_utilities.GetRankZeroDataCommunicator())

            KratosMappingDataTransferOperator.__dummy_model = model
            KratosMappingDataTransferOperator.__rank_zero_model_part = rank_zero_model_part

        return KratosMappingDataTransferOperator.__rank_zero_model_part
