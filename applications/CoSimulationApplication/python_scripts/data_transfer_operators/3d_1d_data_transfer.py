# Importing the base class
from KratosMultiphysics.CoSimulationApplication.base_classes.co_simulation_data_transfer_operator import CoSimulationDataTransferOperator

# Importing the Kratos Library
import KratosMultiphysics as KM
import KratosMultiphysics.MappingApplication as KratosMapping

# CoSimulation imports
import KratosMultiphysics.CoSimulationApplication as KratosCoSim
import KratosMultiphysics.CoSimulationApplication.co_simulation_tools as cs_tools
import KratosMultiphysics.CoSimulationApplication.colors as colors

def Create(*args):
    return Kratos3D1DDataTransferOperator(*args)

class Kratos3D1DDataTransferOperator(CoSimulationDataTransferOperator):
    """DataTransferOperator that transfers values from one 3D interface (ModelPart) to another 1D interface (ModelPart).
    The data_transfer_3d_1ds of the Kratos-MappingApplication are used
    """

    def __init__(self, settings, parent_coupled_solver_data_communicator):
        """
        Initializes the data transfer object with the given settings and sets up
        internal structures for managing data transfer processes between solvers.

        Parameters:
        ----------
        settings : KratosMultiphysics.Parameters
            A Parameters object containing the configuration for data transfer. 
            Must include `"3d_1d_data_transfer_settings"` to specify relevant 
            transfer settings.

        parent_coupled_solver_data_communicator : CoupledSolverDataCommunicator
            The communicator object that facilitates data sharing between solvers 
            in the coupled analysis.

        Raises:
        ------
        Exception
            If the function is not called within an MPI distributed environment.
        Exception
            If the required `"3d_1d_data_transfer_settings"` is missing in `settings`.

        Notes:
        ------
        - This constructor is intended for usage in an MPI environment only. It 
        checks for distributed mode (`KM.IsDistributedRun()`), ensuring the 
        function is not misused outside of MPI.
        - The `__data_transfer_process` dictionary is initialized as an empty 
        structure for storing and managing the data transfer processes.
        - The `origin_is_3d` attribute is set to `None` initially and will be determined 
        based on the geometry of the model part in subsequent processes.
        """
        if KM.IsDistributedRun():
            raise Exception("This function can only be called when Kratos is running in MPI!")

        if not settings.Has("3d_1d_data_transfer_settings"):
            raise Exception('No "3d_1d_data_transfer_settings" provided!')

        # Initialize the parent class with the provided settings and communicator
        super().__init__(settings, parent_coupled_solver_data_communicator)

        # Initialize attributes for tracking geometry dimension and data transfer processes
        self.origin_is_3d = None
        self.__data_transfer_process = {}

    def _ExecuteTransferData(self, from_solver_data, to_solver_data, transfer_options):
        """
        Executes the data transfer between two solvers, either by retrieving and using an existing 
        data transfer process or creating a new one if none exists. The function also allows for 
        specific transfer options such as toggling the sign of the transfer data.

        Parameters:
        ----------
        from_solver_data : SolverData
            Data object containing information about the origin solver's model part and variable to transfer.

        to_solver_data : SolverData
            Data object containing information about the destination solver's model part and variable to receive data.

        transfer_options : TransferOptions
            Additional options for data transfer, such as swapping the sign of transferred values.

        Notes:
        ------
        The function uses internal dictionaries to manage transfer processes, utilizing unique identifiers 
        based on model part names and variables to ensure consistent data mapping between solvers.

        Steps:
        ------
        1. **Identify or Create Transfer Process**:
        - Identifiers for the origin and destination model parts are generated, combining solver names, 
            model part names, and variable names.
        - If a transfer process already exists between the identified model parts, it is reused.
        - If no transfer process exists, a new one is created with the provided settings and options.

        2. **Set Transfer Parameters**:
        - Custom settings are created and added to manage specific transfer configurations, 
            such as swapping signs or adjusting for 3D/1D transfers.

        3. **Execute Transfer**:
        - Depending on whether the origin is in 3D, specific transpose settings are applied, and the 
            transfer is executed.

        Raises:
        ------
        ValueError
            If the `transfer_options` provided are invalid or unsupported.
        """
        model_part_origin_name = from_solver_data.model_part_name
        variable_origin        = from_solver_data.variable
        identifier_origin      = from_solver_data.solver_name + "." + model_part_origin_name + "." + variable_origin.Name()

        model_part_destination_name = to_solver_data.model_part_name
        variable_destination        = to_solver_data.variable
        identifier_destination      = to_solver_data.solver_name + "." + model_part_destination_name + "." + variable_destination.Name()

        identifier_tuple         = (identifier_origin, identifier_destination)
        inverse_identifier_tuple = (identifier_destination, identifier_origin)

        if identifier_tuple in self.__data_transfer_process:
            if self.origin_is_3d == True:
                self.__data_transfer_process[identifier_tuple].Set(KM.Mapper.USE_TRANSPOSE)
            self.__data_transfer_process[identifier_tuple].Execute()
            if self.origin_is_3d == True:
                self.__data_transfer_process[identifier_tuple].Reset(KM.Mapper.USE_TRANSPOSE)
        elif inverse_identifier_tuple in self.__data_transfer_process:
            if self.origin_is_3d == False:
                self.__data_transfer_process[inverse_identifier_tuple].Set(KM.Mapper.USE_TRANSPOSE)
            self.__data_transfer_process[inverse_identifier_tuple].Execute()
            if self.origin_is_3d == False:
                self.__data_transfer_process[inverse_identifier_tuple].Reset(KM.Mapper.USE_TRANSPOSE)
        else:
            model_part_origin      = from_solver_data.GetModelPart()
            model_part_destination = to_solver_data.GetModelPart()

            # Logging information if echo level is set
            if self.echo_level > 0:
                info_msg = (
                    f"Creating Mapper:\n    Origin: ModelPart '{model_part_origin_name}' "
                    f"of solver '{from_solver_data.solver_name}'\n"
                    f"    Destination: ModelPart '{model_part_destination_name}' "
                    f"of solver '{to_solver_data.solver_name}'"
                )
                cs_tools.cs_print_info(colors.bold(self._ClassName()), info_msg)

            # Creating parameters for the data transfer
            # TODO: I should check if the parameters change between different data transfers
            parameters = KM.Parameters(self.settings["3d_1d_data_transfer_settings"].WriteJsonString())
            if not parameters.Has("origin_variables"):
                parameters.AddEmptyArray("origin_variables")
            parameters["origin_variables"].Append(variable_origin.Name())
            if not parameters.Has("destination_variables"):
                parameters.AddEmptyArray("destination_variables")
            parameters["destination_variables"].Append(variable_destination.Name())
            for transfer_option in transfer_options.GetStringArray():
                if transfer_option == "swap_sign":
                    parameters["swap_sign"].SetBool(True)

            # Check if the origin model part is 3D or 1D
            self.origin_is_3d = self.__check_model_part_3D(model_part_origin)
            # Clone is necessary because the settings are validated and defaults assigned, which could influence the creation of other data transfers
            self.__data_transfer_process[identifier_tuple] = KratosCoSim.DataTransfer3D1DProcess(model_part_origin, model_part_destination, parameters.Clone())
            self.__data_transfer_process[inverse_identifier_tuple] = KratosCoSim.DataTransfer3D1DProcess(model_part_origin, model_part_destination, parameters.Clone())

            # Execute the data transfer
            if self.origin_is_3d:
                self.__data_transfer_process[identifier_tuple].Set(KM.Mapper.USE_TRANSPOSE)
            self.__data_transfer_process[identifier_tuple].Execute()
            if self.origin_is_3d:
                self.__data_transfer_process[identifier_tuple].Reset(KM.Mapper.USE_TRANSPOSE)

    def _Check(self, from_solver_data, to_solver_data):
        """
        Validates that the `from_solver_data` and `to_solver_data` contain only 
        nodal data for transfer operations. Raises an exception if the data 
        location is not at the nodes.

        Parameters:
        ----------
        from_solver_data : SolverData
            The source data to check, containing information about the solver and 
            data location in the model part.

        to_solver_data : SolverData
            The destination data to check, containing information about the solver 
            and data location in the model part.

        Raises:
        ------
        Exception
            If the location in `from_solver_data` or `to_solver_data` is not at 
            the nodes, an exception is raised specifying the issue, the model part, 
            and the solver name.

        Notes:
        ------
        This function assumes `from_solver_data` and `to_solver_data` are instances 
        containing a `location` attribute to specify where the data is located, as 
        well as `model_part_name` and `solver_name` attributes for error reporting.
        """
        def CheckData(data_to_check):
            if "node" not in data_to_check.location:
                raise Exception(
                    f'Transfer only supports nodal values! "{self._ClassName()}"\n'
                    f'Checking ModelPart "{data_to_check.model_part_name}" of solver "{data_to_check.solver_name}"'
                )

        # Perform the nodal data location check on both the source and destination data
        CheckData(from_solver_data)
        CheckData(to_solver_data)

    def __check_model_part_3D(self, model_part):
        """
        Checks whether the given `model_part` consists entirely of 3D elements and 
        verifies that all elements have the same geometry type. 

        Parameters:
        ----------
        model_part : ModelPart
            The model part to check, which should contain elements representing a 3D structure.
        Returns:
        -------
        bool
            True if all elements are 3D and have the same geometry type, False if any 
            element is 1D.
        Raises:
        ------
        ValueError
            If elements in `model_part` have different geometry types, an exception is raised.
        """
        first_geometry_type = None

        for elem in model_part.Elements:
            geom = elem.GetGeometry()
            # Check for consistent geometry type
            if first_geometry_type is None:
                first_geometry_type = geom.LocalSpaceDimension()
            elif geom.LocalSpaceDimension() != first_geometry_type:
                raise ValueError("Inconsistent geometry types found in model_part.")

        # Check if the element is 1D
        if first_geometry_type == 1:
            return False
        else:
            return True

    @classmethod
    def _GetDefaultParameters(cls):
        this_defaults = KM.Parameters("""{
            "3d_1d_data_transfer_settings" : {
                "origin_variables"         : [],
                "destination_variables"    : [],
                "swap_sign"                : false,
                "interpolate_parameters"   : {
                    "data_transfer_3d_1d_type" : "nearest_element",
                    "echo_level"  : 0,
                    "search_settings" : {
                        "max_num_search_iterations"     : 8,
                        "echo_level"                    : 0
                    }
                },
                "search_parameters"        :  {
                    "allocation_size"         : 100,
                    "bucket_size"             : 4,
                    "search_factor"           : 2.0,
                    "search_increment_factor" : 1.5
                }
            }
        }""")
        this_defaults.AddMissingParameters(super()._GetDefaultParameters())
        return this_defaults

    @classmethod
    def _GetListAvailableTransferOptions(cls):
        return ["swap_sign"]
