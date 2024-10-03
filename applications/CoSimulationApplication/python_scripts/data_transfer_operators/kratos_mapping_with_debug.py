# Importing the base class
from KratosMultiphysics.CoSimulationApplication.data_transfer_operators.kratos_mapping import KratosMappingDataTransferOperator

# Importing the Kratos Library
import KratosMultiphysics as KM
from KratosMultiphysics.MappingApplication import python_mapper_factory

# CoSimulation imports
import KratosMultiphysics.CoSimulationApplication.co_simulation_tools as cs_tools
import KratosMultiphysics.CoSimulationApplication.colors as colors

# Import the VTK output process module from KratosMultiphysics
import KratosMultiphysics.vtk_output_process as vtk_output_process

# Import the time module
from time import time

# Import os module
import os

# Import the copy module
import copy

def Create(*args):
    return KratosMappingDataTransferOperatorWithDebug(*args)

class KratosMappingDataTransferOperatorWithDebug(KratosMappingDataTransferOperator):
    """DataTransferOperator that maps values from one interface (ModelPart) to another with VTK output for debugging.
    The mappers of the Kratos-MappingApplication are used
    """
    def __init__(self, settings, parent_coupled_solver_data_communicator):
        # Set the default value for the auxiliary debug counter
        self.auxiliary_debug_counter = False
        if settings.Has("debug_settings"):
            if settings["debug_settings"].Has("auxiliary_debug_counter"):
                self.auxiliary_debug_counter = settings["debug_settings"]["auxiliary_debug_counter"].GetBool()
            # Remove the debug settings from the settings
            settings.RemoveValue("debug_settings")

        # Call the base class constructor
        super().__init__(settings, parent_coupled_solver_data_communicator)

        # Initialize the debug VTK values
        self.counter = 0
        self.sub_counter = 0
        self.number_of_mappers = 0

    def _ExecuteTransferData(self, from_solver_data, to_solver_data, transfer_options):
        # Prepare the solver data
        solver_data = self._PrepareSolverData(from_solver_data, to_solver_data, transfer_options)

        if solver_data.identifier_tuple in self._mappers:
            # Generate VTK output for debugging
            if self.auxiliary_debug_counter:
                self.__GenerateProcessVTK(from_solver_data, to_solver_data, transfer_options)
            self.__PostProcessVTKPre(solver_data.identifier_tuple, from_solver_data, to_solver_data, transfer_options)
            self._mappers[solver_data.identifier_tuple].Map(solver_data.variable_origin, solver_data.variable_destination, solver_data.mapper_flags)
            self.__PostProcessVTKPost(solver_data.identifier_tuple, from_solver_data, to_solver_data, transfer_options)
        elif solver_data.inverse_identifier_tuple in self._mappers:
            # Generate VTK output for debugging
            if self.auxiliary_debug_counter:
                self.__GenerateProcessVTK(from_solver_data, to_solver_data, transfer_options)
            self.__PostProcessVTKPre(solver_data.inverse_identifier_tuple, from_solver_data, to_solver_data, transfer_options)
            self._mappers[solver_data.inverse_identifier_tuple].InverseMap(solver_data.variable_destination, solver_data.variable_origin, solver_data.mapper_flags)
            self.__PostProcessVTKPost(solver_data.inverse_identifier_tuple, from_solver_data, to_solver_data, transfer_options)
        else:
            # Generate VTK output for debugging
            self.__GenerateProcessVTK(from_solver_data, to_solver_data, transfer_options)

            # Get the model parts
            model_part_origin      = self._GetModelPartFromInterfaceData(from_solver_data)
            model_part_destination = self._GetModelPartFromInterfaceData(to_solver_data)

            if model_part_origin.IsDistributed() or model_part_destination.IsDistributed():
                mapper_create_fct = python_mapper_factory.CreateMPIMapper
            else:
                mapper_create_fct = python_mapper_factory.CreateMapper

            if self.echo_level > 0:
                info_msg  = "Creating Mapper:\n"
                info_msg += '    Origin: ModelPart "{}" of solver "{}"\n'.format(solver_data.model_part_origin_name, from_solver_data.solver_name)
                info_msg += '    Destination: ModelPart "{}" of solver "{}"'.format(solver_data.model_part_destination_name, to_solver_data.solver_name)

                cs_tools.cs_print_info(colors.bold(self._ClassName()), info_msg)

            mapper_creation_start_time = time()
            self._mappers[solver_data.identifier_tuple] = mapper_create_fct(model_part_origin, model_part_destination, self.settings["mapper_settings"].Clone()) # Clone is necessary because the settings are validated and defaults assigned, which could influence the creation of other mappers
            self.number_of_mappers += 1

            if self.echo_level > 2:
                cs_tools.cs_print_info(colors.bold(self._ClassName()), "Creating Mapper took: {0:.{1}f} [s]".format(time()-mapper_creation_start_time,2))

            self.__PostProcessVTKPre(solver_data.identifier_tuple, from_solver_data, to_solver_data, transfer_options)
            self._mappers[solver_data.identifier_tuple].Map(solver_data.variable_origin, solver_data.variable_destination, solver_data.mapper_flags)
            self.__PostProcessVTKPost(solver_data.identifier_tuple, from_solver_data, to_solver_data, transfer_options)

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
                "output_sub_model_parts"             : true,
                "nodal_solution_step_data_variables" : [],
                "nodal_data_value_variables"         : [],
                "output_interval"                    : 1,
                "output_path"                        : ""
            }
        }""")
        return vtk_default_output_parameters

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
        model_part_origin, model_part_destination, model_part_origin_name, model_part_destination_name, variable_origin, variable_destination, mapper_flags, identifier_origin, identifier_destination, identifier_tuple, inverse_identifier_tuple = self._PrepareSolverData(from_solver_data, to_solver_data, transfer_options)

        # Get the model associated with the origin model part
        model_origin = model_part_origin.GetModel()
        model_destination = model_part_destination.GetModel()

        # Create VTK output settings for the origin model part
        origin_settings = self.__DefaultVTKSettings()
        origin_settings["Parameters"]["model_part_name"].SetString(model_part_origin_name)
        origin_settings["Parameters"]["output_path"].SetString("debug_vtk_output_" + variable_origin.Name() + "_" + identifier_origin)
        if mapper_flags.Is(KM.Mapper.FROM_NON_HISTORICAL):
            origin_settings["Parameters"]["nodal_data_value_variables"].Append(variable_origin.Name())
        else:
            origin_settings["Parameters"]["nodal_solution_step_data_variables"].Append(variable_origin.Name())

        # Create VTK output settings for the destination model part
        destination_settings = self.__DefaultVTKSettings()
        destination_settings["Parameters"]["model_part_name"].SetString(model_part_destination_name)
        destination_settings["Parameters"]["output_path"].SetString("debug_vtk_output_" + variable_destination.Name() + "_" + identifier_destination)
        if mapper_flags.Is(KM.Mapper.TO_NON_HISTORICAL):
            destination_settings["Parameters"]["nodal_data_value_variables"].Append(variable_destination.Name())
        else:
            destination_settings["Parameters"]["nodal_solution_step_data_variables"].Append(variable_destination.Name())

        return model_origin, model_destination, origin_settings, destination_settings, variable_origin.Name(), variable_destination.Name(), identifier_tuple, inverse_identifier_tuple

    def __GenerateProcessVTK(self, from_solver_data, to_solver_data, transfer_options):
        """
        Execute VTK output processing using provided settings on the specified model part.

        Args:
            self: The instance of the class.
            from_solver_data (CoSimulationData): The data from the solver to transfer from.
            to_solver_data (CoSimulationData): The data from the solver to transfer to.
            transfer_options (KM.Flags): The flags to be used for the mapping.
        """
        # Define the name appendix for the VTK output
        name_prefix = ""

        # Increase the counter for debugging
        if self.auxiliary_debug_counter:
            self.sub_counter += 1
            if self.sub_counter >= 2 * self.number_of_mappers: # 2 because of pre and post mapping
                self.sub_counter = 0
                self.counter += 1
            name_prefix = f"MAPPING_{self.counter}"

        # Prepare the settings for VTK output processing
        model_origin, model_destination, origin_settings, destination_settings, variable_origin_name, variable_destination_name, identifier_tuple, inverse_identifier_tuple = self.__PrepareSettings(from_solver_data, to_solver_data, transfer_options)
        variable_identifier_tuple = (variable_origin_name, variable_destination_name)

        # Define __debug_vtk_pre if not defined
        if not hasattr(self, "_KratosMappingDataTransferOperatorWithDebug__debug_vtk_pre"):
            self.__debug_vtk_pre = {}

        # Define __debug_vtk_post if not defined
        if not hasattr(self, "_KratosMappingDataTransferOperatorWithDebug__debug_vtk_post"):
            self.__debug_vtk_post = {}

        # Check if the VTK output processes are already defined
        current_identifier = None
        if identifier_tuple in self.__debug_vtk_pre:
            current_identifier = identifier_tuple
        elif inverse_identifier_tuple in self.__debug_vtk_pre:
            current_identifier = inverse_identifier_tuple
        if current_identifier is None:
            self.__debug_vtk_pre[identifier_tuple] = {}
            self.__debug_vtk_post[identifier_tuple] = {}
            current_identifier = identifier_tuple

        # Get the output paths
        origin_output_path = origin_settings["Parameters"]["output_path"].GetString()
        destination_output_path = destination_settings["Parameters"]["output_path"].GetString()

        # Set the output paths
        pre_origin_settings = copy.deepcopy(origin_settings)
        path_splitted = os.path.split(origin_output_path)
        directory_name = path_splitted[0]
        for i in range(1, len(path_splitted)-1):
            directory_name = os.path.join(directory_name, path_splitted[i])
        base_name = path_splitted[-1]
        if self.auxiliary_debug_counter:
            directory_name = os.path.join(directory_name, name_prefix)
        # Ensure the directory exists
        if directory_name != "":
            os.makedirs(directory_name, exist_ok=True)
        file_path = os.path.join(directory_name, "pre_map_" + base_name)
        pre_origin_settings["Parameters"]["output_path"].SetString(file_path)

        # Create a VTK output process object with the provided settings and the current model
        pre_process_origin = vtk_output_process.Factory(pre_origin_settings, model_origin)
        # Initialize and execute the VTK output process
        pre_process_origin.ExecuteInitialize()
        pre_process_origin.ExecuteInitializeSolutionStep()

        # Set the output paths
        pre_destination_settings = copy.deepcopy(destination_settings)
        path_splitted = os.path.split(destination_output_path)
        directory_name = path_splitted[0]
        for i in range(1, len(path_splitted)-1):
            directory_name = os.path.join(directory_name, path_splitted[i])
        base_name = path_splitted[-1]
        if self.auxiliary_debug_counter:
            directory_name = os.path.join(directory_name, name_prefix)
        file_path = os.path.join(directory_name, "pre_map_" + base_name)
        pre_destination_settings["Parameters"]["output_path"].SetString(file_path)

        # Create a VTK output process object with the provided settings and the current model
        pre_process_destination = vtk_output_process.Factory(pre_destination_settings, model_destination)
        # Initialize and execute the VTK output process
        pre_process_destination.ExecuteInitialize()
        pre_process_destination.ExecuteInitializeSolutionStep()

        # Store the VTK output processes for debugging
        self.__debug_vtk_pre[current_identifier][variable_identifier_tuple] = (pre_process_origin, pre_process_destination)

        # Set the output paths
        post_origin_settings = copy.deepcopy(origin_settings)
        path_splitted = os.path.split(origin_output_path)
        directory_name = path_splitted[0]
        for i in range(1, len(path_splitted)-1):
            directory_name = os.path.join(directory_name, path_splitted[i])
        base_name = path_splitted[-1]
        if self.auxiliary_debug_counter:
            directory_name = os.path.join(directory_name, name_prefix)
        file_path = os.path.join(directory_name, "post_map_" + base_name)
        post_origin_settings["Parameters"]["output_path"].SetString(file_path)

        # Create a VTK output process object with the provided settings and the current model
        post_process_origin = vtk_output_process.Factory(post_origin_settings, model_origin)
        # Initialize and execute the VTK output process
        post_process_origin.ExecuteInitialize()
        post_process_origin.ExecuteInitializeSolutionStep()

        # Set the output paths
        post_destination_settings = copy.deepcopy(destination_settings)
        path_splitted = os.path.split(destination_output_path)
        directory_name = path_splitted[0]
        for i in range(1, len(path_splitted)-1):
            directory_name = os.path.join(directory_name, path_splitted[i])
        base_name = path_splitted[-1]
        if self.auxiliary_debug_counter:
            directory_name = os.path.join(directory_name, name_prefix)
        file_path = os.path.join(directory_name, "post_map_" + base_name)
        post_destination_settings["Parameters"]["output_path"].SetString(file_path)

        # Create a VTK output process object with the provided settings and the current model
        post_process_destination = vtk_output_process.Factory(post_destination_settings, model_destination)
        # Initialize and execute the VTK output process
        post_process_destination.ExecuteInitialize()
        post_process_destination.ExecuteInitializeSolutionStep()

        # Store the VTK output processes for debugging
        self.__debug_vtk_post[current_identifier][variable_identifier_tuple] = (post_process_origin, post_process_destination)

    def __GetIdentifier(self, identifier_tuple, from_solver_data, to_solver_data, transfer_options):
        """
        Gets the identifier used for the postprocess.

        Args:
            self: The instance of the class.
            identifier_tuple (tuple): Tuple containing the identifiers of the model parts to be processed.
            from_solver_data (CoSimulationData): The data from the solver to transfer from.
            to_solver_data (CoSimulationData): The data from the solver to transfer to.
            transfer_options (KM.Flags): The flags to be used for the mapping.
        """
        variable_origin = from_solver_data.variable
        variable_destination = to_solver_data.variable
        variable_identifier_tuple = (variable_origin.Name(), variable_destination.Name())

        inverse_identifier_tuple = (identifier_tuple[1], identifier_tuple[0])
        if identifier_tuple in self.__debug_vtk_pre:
            current_identifier = identifier_tuple
        elif inverse_identifier_tuple in self.__debug_vtk_pre:
            current_identifier = inverse_identifier_tuple
        else:
            # Generate VTK output for debugging
            self.__GenerateProcessVTK(from_solver_data, to_solver_data, transfer_options)
            if identifier_tuple in self.__debug_vtk_pre:
                current_identifier = identifier_tuple
            elif inverse_identifier_tuple in self.__debug_vtk_pre:
                current_identifier = inverse_identifier_tuple
            else:
                raise Exception("Could not generate VTK output for debugging!")

        # Generate VTK output for debugging if not already done
        if not variable_identifier_tuple in self.__debug_vtk_pre[current_identifier]:
            # Generate VTK output for debugging
            self.__GenerateProcessVTK(from_solver_data, to_solver_data, transfer_options)

        return current_identifier, variable_identifier_tuple

    def __PostProcessVTKPre(self, identifier_tuple, from_solver_data, to_solver_data, transfer_options):
        """
        Execute VTK output processing using provided settings on the specified model part. Pre mapping.

        Args:
            self: The instance of the class.
            identifier_tuple (tuple): Tuple containing the identifiers of the model parts to be processed.
            from_solver_data (CoSimulationData): The data from the solver to transfer from.
            to_solver_data (CoSimulationData): The data from the solver to transfer to.
            transfer_options (KM.Flags): The flags to be used for the mapping.
        """
        current_identifier, variable_identifier_tuple = self.__GetIdentifier(identifier_tuple, from_solver_data, to_solver_data, transfer_options)

        self.__debug_vtk_pre[current_identifier][variable_identifier_tuple][0].ExecuteFinalizeSolutionStep()
        self.__debug_vtk_pre[current_identifier][variable_identifier_tuple][0].PrintOutput()
        self.__debug_vtk_pre[current_identifier][variable_identifier_tuple][1].ExecuteFinalizeSolutionStep()
        self.__debug_vtk_pre[current_identifier][variable_identifier_tuple][1].PrintOutput()

    def __PostProcessVTKPost(self, identifier_tuple, from_solver_data, to_solver_data, transfer_options):
        """
        Execute VTK output processing using provided settings on the specified model part. Post mapping.

        Args:
            self: The instance of the class.
            identifier_tuple (tuple): Tuple containing the identifiers of the model parts to be processed.
            from_solver_data (CoSimulationData): The data from the solver to transfer from.
            to_solver_data (CoSimulationData): The data from the solver to transfer to.
            transfer_options (KM.Flags): The flags to be used for the mapping.
        """
        current_identifier, variable_identifier_tuple = self.__GetIdentifier(identifier_tuple, from_solver_data, to_solver_data, transfer_options)

        self.__debug_vtk_post[current_identifier][variable_identifier_tuple][0].ExecuteFinalizeSolutionStep()
        self.__debug_vtk_post[current_identifier][variable_identifier_tuple][0].PrintOutput()
        self.__debug_vtk_post[current_identifier][variable_identifier_tuple][1].ExecuteFinalizeSolutionStep()
        self.__debug_vtk_post[current_identifier][variable_identifier_tuple][1].PrintOutput()