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

# Import the VTK output process module from KratosMultiphysics
import KratosMultiphysics.vtk_output_process as vtk_output_process

# Import the time module
from time import time

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
        self.debug_vtk = self.settings["debug_vtk"].GetBool()
        self.__mappers = {}

    def _ExecuteTransferData(self, from_solver_data, to_solver_data, transfer_options):
        # Prepare the solver data
        model_part_origin, model_part_destination, model_part_origin_name, model_part_destination_name, variable_origin, variable_destination, mapper_flags, identifier_origin, identifier_destination, identifier_tuple, inverse_identifier_tuple = self.__PrepareSolverData(from_solver_data, to_solver_data, transfer_options)

        if identifier_tuple in self.__mappers:
            self.__PostProcessVTK(identifier_tuple, from_solver_data, to_solver_data, transfer_options, True)
            self.__mappers[identifier_tuple].Map(variable_origin, variable_destination, mapper_flags)
            self.__PostProcessVTK(identifier_tuple, from_solver_data, to_solver_data, transfer_options, False)
        elif inverse_identifier_tuple in self.__mappers:
            self.__PostProcessVTK(inverse_identifier_tuple, from_solver_data, to_solver_data, transfer_options, True)
            self.__mappers[inverse_identifier_tuple].InverseMap(variable_destination, variable_origin, mapper_flags)
            self.__PostProcessVTK(inverse_identifier_tuple, from_solver_data, to_solver_data, transfer_options, False)
        else:
            # Get the model parts
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
            self.__PostProcessVTK(identifier_tuple, from_solver_data, to_solver_data, transfer_options, True)
            self.__mappers[identifier_tuple].Map(variable_origin, variable_destination, mapper_flags)
            self.__PostProcessVTK(identifier_tuple, from_solver_data, to_solver_data, transfer_options, False)

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

    def __DefaultVTKSettings(self):
        """
        Create and return a default VTK output settings using KratosMultiphysics (KM) Parameters.

        Args:
            self: The instance of the class.

        Returns:
            KM.Parameters: A parameters object initialized with default settings for VTK output.
        """
        # Define default settings for VTK output
        vtk_default_output_parameters = KM.Parameters("""{
            "Parameters" : {
                "model_part_name"                    : "",
                "file_format"                        : "ascii",
                "entity_type"                        : "element",
                "nodal_solution_step_data_variables" : [],
                "nodal_data_value_variables"         : [],
                "output_interval"                    : 1,
                "output_path"                        : ""
            }
        }""")
        return vtk_default_output_parameters

    def __PrepareSolverData(self, from_solver_data, to_solver_data, transfer_options):
        """
        Prepare the solver data

        Args:
            self: The instance of the class.
            from_solver_data (CoSimulationData): The data from the solver to transfer from.
            to_solver_data (CoSimulationData): The data from the solver to transfer to.
            transfer_options (KM.Flags): The flags to be used for the mapping.
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
        model_part_origin      = self.__GetModelPartFromInterfaceData(from_solver_data)
        model_part_destination = self.__GetModelPartFromInterfaceData(to_solver_data)

        return model_part_origin, model_part_destination, model_part_origin_name, model_part_destination_name, variable_origin, variable_destination, mapper_flags, identifier_origin, identifier_destination, identifier_tuple, inverse_identifier_tuple

    def __PrepareSettings(self, from_solver_data, to_solver_data, transfer_options):
        """
        Prepare the settings for VTK output processing.

        Args:
            self: The instance of the class.
            from_solver_data (CoSimulationData): The data from the solver to transfer from.
            to_solver_data (CoSimulationData): The data from the solver to transfer to.
            transfer_options (KM.Flags): The flags to be used for the mapping.
        """
        # Prepare the solver data
        model_part_origin, model_part_destination, model_part_origin_name, model_part_destination_name, variable_origin, variable_destination, mapper_flags, identifier_origin, identifier_destination, identifier_tuple, inverse_identifier_tuple = self.__PrepareSolverData(from_solver_data, to_solver_data, transfer_options)

        # Get the model associated with the origin model part
        model_origin = model_part_origin.GetModel()
        model_destination = model_part_destination.GetModel()

        # Create VTK output settings for the origin model part
        origin_settings = self.__DefaultVTKSettings()
        origin_settings["Parameters"]["model_part_name"].SetString(model_part_origin_name)
        origin_settings["Parameters"]["output_path"].SetString("debug_vtk_output_" + identifier_origin)
        if mapper_flags.Is(KM.Mapper.FROM_NON_HISTORICAL):
            origin_settings["Parameters"]["nodal_data_value_variables"].Append(variable_origin.Name())
        else:
            origin_settings["Parameters"]["nodal_solution_step_data_variables"].Append(variable_origin.Name())
        if model_part_origin.NumberOfElements() == 0:
            origin_settings["Parameters"]["entity_type"].SetString("condition")

        # Create VTK output settings for the destination model part
        destination_settings = self.__DefaultVTKSettings()
        destination_settings["Parameters"]["model_part_name"].SetString(model_part_destination_name)
        destination_settings["Parameters"]["output_path"].SetString("debug_vtk_output_" + identifier_destination)
        if mapper_flags.Is(KM.Mapper.TO_NON_HISTORICAL):
            destination_settings["Parameters"]["nodal_data_value_variables"].Append(variable_destination.Name())
        else:
            destination_settings["Parameters"]["nodal_solution_step_data_variables"].Append(variable_destination.Name())
        if model_part_destination.NumberOfElements() == 0:
            destination_settings["Parameters"]["entity_type"].SetString("condition")

        return model_origin, model_destination, origin_settings, destination_settings, identifier_tuple

    def __GenerateProcessVTK(self, from_solver_data, to_solver_data, transfer_options):
        """
        Execute VTK output processing using provided settings on the specified model part.

        Args:
            self: The instance of the class.
            from_solver_data (CoSimulationData): The data from the solver to transfer from.
            to_solver_data (CoSimulationData): The data from the solver to transfer to.
            transfer_options (KM.Flags): The flags to be used for the mapping.
        """
        if self.debug_vtk:
            model_origin, model_destination, origin_settings, destination_settings, identifier_tuple = self.__PrepareSettings(from_solver_data, to_solver_data, transfer_options)

            # Get the output paths
            origin_output_path = origin_settings["Parameters"]["output_path"].GetString()
            destination_output_path = destination_settings["Parameters"]["output_path"].GetString()

            # Copy settings for pre and post mapping
            pre_origin_settings = KM.Parameters(origin_settings.PrettyPrintJsonString())
            post_origin_settings = KM.Parameters(origin_settings.PrettyPrintJsonString())
            pre_destination_settings = KM.Parameters(destination_settings.PrettyPrintJsonString())
            post_destination_settings = KM.Parameters(destination_settings.PrettyPrintJsonString())

            # Define __debug_vtk_settings if not defined
            if not hasattr(self, "__debug_vtk_settings"):
                self.__debug_vtk_settings = {}
            if not hasattr(self, "__debug_vtk_variables"):
                self.__debug_vtk_variables = {}
            self.__debug_vtk_settings[identifier_tuple] = (pre_origin_settings, post_origin_settings, pre_destination_settings, post_destination_settings)
            self.__debug_vtk_variables[identifier_tuple] = [(from_solver_data.variable, to_solver_data.variable)]

            # Set the output paths
            pre_origin_settings["Parameters"]["output_path"].SetString(origin_output_path + "_pre_map")

            # Create a VTK output process object with the provided settings and the current model
            pre_process_origin = vtk_output_process.Factory(pre_origin_settings, model_origin)
            # Initialize and execute the VTK output process
            pre_process_origin.ExecuteInitialize()
            pre_process_origin.ExecuteInitializeSolutionStep()

            # Set the output paths
            post_origin_settings["Parameters"]["output_path"].SetString(origin_output_path + "_post_map")

            # Create a VTK output process object with the provided settings and the current model
            post_process_origin = vtk_output_process.Factory(post_origin_settings, model_origin)
            # Initialize and execute the VTK output process
            post_process_origin.ExecuteInitialize()
            post_process_origin.ExecuteInitializeSolutionStep()

            # Set the output paths
            pre_destination_settings["Parameters"]["output_path"].SetString(destination_output_path + "_pre_map")

            # Create a VTK output process object with the provided settings and the current model
            pre_process_destination = vtk_output_process.Factory(pre_destination_settings, model_destination)
            # Initialize and execute the VTK output process
            pre_process_destination.ExecuteInitialize()
            pre_process_destination.ExecuteInitializeSolutionStep()

            # Set the output paths
            post_destination_settings["Parameters"]["output_path"].SetString(destination_output_path + "_post_map")

            # Create a VTK output process object with the provided settings and the current model
            post_process_destination = vtk_output_process.Factory(post_destination_settings, model_destination)
            # Initialize and execute the VTK output process
            post_process_destination.ExecuteInitialize()
            post_process_destination.ExecuteInitializeSolutionStep()

            # Define __debug_vtk if not defined
            if not hasattr(self, "__debug_vtk"):
                self.__debug_vtk = {}

            # Store the VTK output processes for debugging
            self.__debug_vtk[identifier_tuple] = (pre_process_origin, post_process_origin, pre_process_destination, post_process_destination)

    def __PostProcessVTK(self, identifier_tuple, from_solver_data, to_solver_data, transfer_options, pre_map = False):
        """
        Execute VTK output processing using provided settings on the specified model part.

        Args:
            self: The instance of the class.
            identifier_tuple (tuple): Tuple containing the identifiers of the model parts to be processed.
            from_solver_data (CoSimulationData): The data from the solver to transfer from.
            to_solver_data (CoSimulationData): The data from the solver to transfer to.
            transfer_options (KM.Flags): The flags to be used for the mapping.
            pre_map (bool): Flag to indicate if the processing is done before or after mapping.
        """
        if self.debug_vtk:
            # Define __debug_vtk if not defined
            if not hasattr(self, "__debug_vtk"):
                self.__debug_vtk = {}
            # Generate VTK output for debugging if not already done
            if not identifier_tuple in self.__debug_vtk:
                # Generate VTK output for debugging
                self.__GenerateProcessVTK(from_solver_data, to_solver_data, transfer_options)

            # Get the current post process identifier
            if identifier_tuple in self.__debug_vtk_variables:
                current_post_process_identifier = identifier_tuple
            elif (identifier_tuple[1], identifier_tuple[0]) in self.__debug_vtk_variables:
                current_post_process_identifier = (identifier_tuple[1], identifier_tuple[0])
            else:
                raise Exception('Definition not found in the VTK output processes!')

            # Prepare the solver data
            model_part_origin, model_part_destination, model_part_origin_name, model_part_destination_name, variable_origin, variable_destination, mapper_flags, identifier_origin, identifier_destination, identifier_tuple, inverse_identifier_tuple = self.__PrepareSolverData(from_solver_data, to_solver_data, transfer_options)

            if (variable_origin, variable_destination) in self.__debug_vtk_variables[current_post_process_identifier]:
                if pre_map:
                    self.__debug_vtk[current_post_process_identifier][0].ExecuteFinalizeSolutionStep()
                    self.__debug_vtk[current_post_process_identifier][0].PrintOutput()
                    self.__debug_vtk[current_post_process_identifier][2].ExecuteFinalizeSolutionStep()
                    self.__debug_vtk[current_post_process_identifier][2].PrintOutput()
                else:
                    self.__debug_vtk[current_post_process_identifier][1].ExecuteFinalizeSolutionStep()
                    self.__debug_vtk[current_post_process_identifier][1].PrintOutput()
                    self.__debug_vtk[current_post_process_identifier][3].ExecuteFinalizeSolutionStep()
                    self.__debug_vtk[current_post_process_identifier][3].PrintOutput()
            else:
                self.__debug_vtk_variables[current_post_process_identifier].append((variable_origin, variable_destination))
                settings = self.__debug_vtk_settings[current_post_process_identifier]
                model_origin = model_part_origin.GetModel()
                model_destination = model_part_destination.GetModel()

                # Update the settings
                pre_origin_settings = settings[0]
                post_origin_settings = settings[2]
                pre_destination_settings = settings[1]
                post_destination_settings = settings[3]

                # Create a VTK output process object with the provided settings and the current model
                model_part_name = pre_origin_settings["Parameters"]["model_part_name"].GetString()
                if model_origin.HasModelPart(model_part_name):
                    if mapper_flags.Is(KM.Mapper.FROM_NON_HISTORICAL):
                        pre_origin_settings["Parameters"]["nodal_data_value_variables"].Append(variable_origin.Name())
                    else:
                        pre_origin_settings["Parameters"]["nodal_solution_step_data_variables"].Append(variable_origin.Name())
                    pre_process_origin = vtk_output_process.Factory(pre_origin_settings, model_origin)
                elif model_destination.HasModelPart(model_part_name):
                    if mapper_flags.Is(KM.Mapper.FROM_NON_HISTORICAL):
                        pre_origin_settings["Parameters"]["nodal_data_value_variables"].Append(variable_destination.Name())
                    else:
                        pre_origin_settings["Parameters"]["nodal_solution_step_data_variables"].Append(variable_destination.Name())
                    pre_process_origin = vtk_output_process.Factory(pre_origin_settings, model_destination)
                else:
                    raise Exception('ModelPart "{}" not found in the model!'.format(model_part_name))
                # Initialize and execute the VTK output process
                pre_process_origin.ExecuteInitialize()
                pre_process_origin.ExecuteInitializeSolutionStep()
                if pre_map:
                    pre_process_origin.ExecuteFinalizeSolutionStep()
                    pre_process_origin.PrintOutput()

                # Create a VTK output process object with the provided settings and the current model
                model_part_name = post_origin_settings["Parameters"]["model_part_name"].GetString()
                if model_origin.HasModelPart(model_part_name):
                    if mapper_flags.Is(KM.Mapper.FROM_NON_HISTORICAL):
                        post_origin_settings["Parameters"]["nodal_data_value_variables"].Append(variable_origin.Name())
                    else:
                        post_origin_settings["Parameters"]["nodal_solution_step_data_variables"].Append(variable_origin.Name())
                    post_process_origin = vtk_output_process.Factory(post_origin_settings, model_origin)
                elif model_destination.HasModelPart(model_part_name):
                    if mapper_flags.Is(KM.Mapper.FROM_NON_HISTORICAL):
                        post_origin_settings["Parameters"]["nodal_data_value_variables"].Append(variable_destination.Name())
                    else:
                        post_origin_settings["Parameters"]["nodal_solution_step_data_variables"].Append(variable_destination.Name())
                    post_process_origin = vtk_output_process.Factory(post_origin_settings, model_destination)
                else:
                    raise Exception('ModelPart "{}" not found in the model!'.format(model_part_name))
                # Initialize and execute the VTK output process
                post_process_origin.ExecuteInitialize()
                post_process_origin.ExecuteInitializeSolutionStep()
                if not pre_map:
                    post_process_origin.ExecuteFinalizeSolutionStep()
                    post_process_origin.PrintOutput()

                # Create a VTK output process object with the provided settings and the current model
                model_part_name = pre_destination_settings["Parameters"]["model_part_name"].GetString()
                if model_origin.HasModelPart(model_part_name):
                    if mapper_flags.Is(KM.Mapper.TO_NON_HISTORICAL):
                        pre_destination_settings["Parameters"]["nodal_data_value_variables"].Append(variable_origin.Name())
                    else:
                        pre_destination_settings["Parameters"]["nodal_solution_step_data_variables"].Append(variable_origin.Name())
                    pre_process_destination = vtk_output_process.Factory(pre_destination_settings, model_origin)
                elif model_destination.HasModelPart(model_part_name):
                    if mapper_flags.Is(KM.Mapper.TO_NON_HISTORICAL):
                        pre_destination_settings["Parameters"]["nodal_data_value_variables"].Append(variable_destination.Name())
                    else:
                        pre_destination_settings["Parameters"]["nodal_solution_step_data_variables"].Append(variable_destination.Name())
                    pre_process_destination = vtk_output_process.Factory(pre_destination_settings, model_destination)
                else:
                    raise Exception('ModelPart "{}" not found in the model!'.format(model_part_name))
                # Initialize and execute the VTK output process
                pre_process_destination.ExecuteInitialize()
                pre_process_destination.ExecuteInitializeSolutionStep()
                if pre_map:
                    pre_process_destination.ExecuteFinalizeSolutionStep()
                    pre_process_destination.PrintOutput()

                # Create a VTK output process object with the provided settings and the current model
                model_part_name = post_destination_settings["Parameters"]["model_part_name"].GetString()
                if model_origin.HasModelPart(model_part_name):
                    if mapper_flags.Is(KM.Mapper.TO_NON_HISTORICAL):
                        post_destination_settings["Parameters"]["nodal_data_value_variables"].Append(variable_origin.Name())
                    else:
                        post_destination_settings["Parameters"]["nodal_solution_step_data_variables"].Append(variable_origin.Name())
                    post_process_destination = vtk_output_process.Factory(post_destination_settings, model_origin)
                elif model_destination.HasModelPart(model_part_name):
                    if mapper_flags.Is(KM.Mapper.TO_NON_HISTORICAL):
                        post_destination_settings["Parameters"]["nodal_data_value_variables"].Append(variable_destination.Name())
                    else:
                        post_destination_settings["Parameters"]["nodal_solution_step_data_variables"].Append(variable_destination.Name())
                    post_process_destination = vtk_output_process.Factory(post_destination_settings, model_destination)
                else:
                    raise Exception('ModelPart "{}" not found in the model!'.format(model_part_name))
                # Initialize and execute the VTK output process
                post_process_destination.ExecuteInitialize()
                post_process_destination.ExecuteInitializeSolutionStep()
                if not pre_map:
                    post_process_destination.ExecuteFinalizeSolutionStep()
                    post_process_destination.PrintOutput()

                # Store the VTK output processes for debugging
                self.__debug_vtk[current_post_process_identifier] = (pre_process_origin, post_process_origin, pre_process_destination, post_process_destination)