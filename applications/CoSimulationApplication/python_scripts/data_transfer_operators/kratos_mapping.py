# Importing the base class
from KratosMultiphysics.CoSimulationApplication.base_classes.co_simulation_data_transfer_operator import CoSimulationDataTransferOperator

# Importing the Kratos Library
import KratosMultiphysics as KM
from KratosMultiphysics.MappingApplication import python_mapper_factory

# CoSimulation imports
import KratosMultiphysics.CoSimulationApplication.co_simulation_tools as cs_tools
import KratosMultiphysics.CoSimulationApplication.colors as colors

# other imports
from KratosMultiphysics.CoSimulationApplication.utilities import data_communicator_utilities
from time import time

import numpy as np

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
        self.__mappers = {}

        ###ADDED BY SEBASTIAN####
        # Initialize snapshot dictionaries for temperature and auxiliary flux
        self.temperature_solid = []
        self.temperature_fluid = []
        self.face_heat_flux_solid = []
        self.face_heat_flux_fluid = []
        origin_var_name = "TEMPERATURE"
        if not KM.KratosGlobals.HasVariable(origin_var_name):
            err_msg = "\'{}\' variable in \'nodal_unknowns\' is not in KratosGlobals. Please check provided value.".format(origin_var_name)
        if not KM.KratosGlobals.GetVariableType(origin_var_name):
            err_msg = "\'{}\' variable in \'nodal_unknowns\' is not double type. Please check provide double type variables (e.g. [\"DISPLACEMENT_X\",\"DISPLACEMENT_Y\"]).".format(origin_var_name)
        self.origin_var_to_store = KM.KratosGlobals.GetVariable(origin_var_name)# TODO: Make it something that the user can give
        destination_var_name = "FACE_HEAT_FLUX"
        if not KM.KratosGlobals.HasVariable(destination_var_name):
            err_msg = "\'{}\' variable in \'nodal_unknowns\' is not in KratosGlobals. Please check provided value.".format(destination_var_name)
        if not KM.KratosGlobals.GetVariableType(destination_var_name):
            err_msg = "\'{}\' variable in \'nodal_unknowns\' is not double type. Please check provide double type variables (e.g. [\"DISPLACEMENT_X\",\"DISPLACEMENT_Y\"]).".format(destination_var_name)
        self.destination_var_to_store = KM.KratosGlobals.GetVariable(destination_var_name)# TODO: Make it something that the user can give

        # Initiliazation flag
        self.solid_init_flag = False
        self.fluid_init_flag = False

        self.counter = 0

    def __ExtractNodeIds(self, model_part):
        node_ids = np.array([node.Id for node in model_part.Nodes])
        return node_ids

    def ExtractAndStoreDataSnapshot(self, from_solver_data, to_solver_data):
        # Obtain the model parts
        model_part_origin = self.__GetModelPartFromInterfaceData(from_solver_data)
        model_part_destination = self.__GetModelPartFromInterfaceData(to_solver_data)

        data_origin = self.__ExtractData(model_part_origin, self.origin_var_to_store)
        data_destination = self.__ExtractData(model_part_destination, self.destination_var_to_store)

        # Store the snapshots
        from_solver_data_model_part_name = model_part_origin.Name
        to_solver_data_model_part_name = model_part_destination.Name
        if from_solver_data.solver_name=="solid" and len(self.temperature_solid)==0:#TODO: this is hardcoded, think of something else
            data_origin_node_ids = self.__ExtractNodeIds(model_part_origin)
            origin_node_ids_name = f"{from_solver_data_model_part_name}_interface_node_ids.npy"
            np.save(origin_node_ids_name, data_origin_node_ids)
            data_destination_node_ids = self.__ExtractNodeIds(model_part_destination)
            destination_node_ids_name = f"{to_solver_data_model_part_name}_interface_node_ids.npy"
            np.save(destination_node_ids_name, data_destination_node_ids)
        if from_solver_data.solver_name=="solid":
            self.temperature_solid.append(data_origin)
            temperature_data_file_name = f"{from_solver_data_model_part_name}_temperature.npy"
            np.save(temperature_data_file_name, np.array(self.temperature_solid))
            self.face_heat_flux_fluid.append(data_destination)
            face_heat_flux_data_file_name = f"{to_solver_data_model_part_name}_face_heat_flux.npy"
            np.save(face_heat_flux_data_file_name, np.array(self.face_heat_flux_fluid))
        elif from_solver_data.solver_name=="fluid":
            self.temperature_fluid.append(data_origin)
            temperature_data_file_name = f"{from_solver_data_model_part_name}_temperature.npy"
            np.save(temperature_data_file_name, np.array(self.temperature_fluid))
            self.face_heat_flux_solid.append(data_destination)
            face_heat_flux_data_file_name = f"{to_solver_data_model_part_name}_face_heat_flux.npy"
            np.save(face_heat_flux_data_file_name, np.array(self.face_heat_flux_solid))

    def __ExtractData(self, model_part, variable):
        """
        Extracts data for a specified variable from a model part.

        :param model_part: The model part from which to extract data.
        :param variable: The variable for which to extract data.
        :return: A NumPy array with the corresponding variable values.
        """
        # Assuming GetSolutionStepValue(variable) for each node returns a scalar or a vector of variable values
        # across time steps. If it's just a single value per time step, this needs to be accumulated in a loop
        # over time steps which is not shown here.

        variable_values = np.array([node.GetSolutionStepValue(variable) for node in model_part.Nodes])

        return variable_values

    def _ExecuteTransferDataWithTransferOperator(self, from_solver_data, to_solver_data, transfer_options):

        self.P_temp_to_temp = np.load("transfer_operator_P_temperature_to_temperature.npy")  # Load the projection operator P temp
        self.P_temp_to_face_heat_flux = np.load("transfer_operator_P_temperature_to_face_heat_flux.npy")  # Load the projection operator P flux

        # Obtain the model parts
        model_part_origin = self.__GetModelPartFromInterfaceData(from_solver_data)
        model_part_destination = self.__GetModelPartFromInterfaceData(to_solver_data)

        if from_solver_data.solver_name=="solid":
            # Retrieve temperatures from the solid interface
            solid_interface_temperatures = np.array([
                node.GetSolutionStepValue(KM.TEMPERATURE)
                for node in model_part_origin.Nodes
            ])

            # Apply the projection operator P to map temperatures
            mapped_fluid_temperatures = np.dot(self.P_temp_to_temp, solid_interface_temperatures)

            # mapped_fluid_temperatures = np.load("GENERIC_Interface_fluid_temperature.npy")[self.counter,:]

            # if not self.solid_init_flag:
            #     mapped_fluid_temperatures *= 0.0
            #     self.solid_init_flag = True

            # Assign the mapped temperatures to the fluid interface nodes
            for i, node in enumerate(model_part_destination.Nodes):
                node.SetSolutionStepValue(KM.TEMPERATURE, mapped_fluid_temperatures[i])

        elif from_solver_data.solver_name=="fluid":
            # Retrieve temperatures from the solid interface
            fluid_interface_temperatures = np.array([
                node.GetSolutionStepValue(KM.TEMPERATURE)
                for node in model_part_origin.Nodes
            ])

            # Apply the projection operator P to map temperatures to face heat flux
            mapped_solid_face_heat_flux = np.dot(self.P_temp_to_face_heat_flux, fluid_interface_temperatures)

            # mapped_solid_face_heat_flux = np.load("GENERIC_Interface_solid_face_heat_flux.npy")[self.counter,:]

            # if not self.fluid_init_flag:
            #     mapped_solid_face_heat_flux *= 0.0
            #     self.fluid_init_flag = True

            # Assign the mapped temperatures to the fluid interface nodes
            for i, node in enumerate(model_part_destination.Nodes):
                node.SetSolutionStepValue(KM.FACE_HEAT_FLUX, mapped_solid_face_heat_flux[i])

        self.counter += 1
    ###ADDED BY SEBASTIAN####

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
        return cls.__mapper_flags_dict.keys()

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
    def __GetModelPartFromInterfaceData(interface_data):
        """If the solver does not exist on this rank, then pass a
        dummy ModelPart to the Mapper that has a DataCommunicator
        that is not defined on this rank
        """
        if interface_data.IsDefinedOnThisRank():
            return interface_data.GetModelPart()
        else:
            return KratosMappingDataTransferOperator.__GetRankZeroModelPart()

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