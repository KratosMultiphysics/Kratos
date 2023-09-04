# Importing the main Kratos Library
import KratosMultiphysics as KM

# Importing base class and utilities
from KratosMultiphysics.CoSimulationApplication.base_classes.co_simulation_solver_wrapper import CoSimulationSolverWrapper
from KratosMultiphysics.CoSimulationApplication.utilities import model_part_utilities
from KratosMultiphysics.CoSimulationApplication.utilities import data_communicator_utilities
import KratosMultiphysics.CoSimulationApplication.co_simulation_tools as cs_tools

# Importing standard required modules
import os
import sys
from pathlib import Path

def Create(settings, model, solver_name):
    """
    Create an instance of the FluxWrapper.

    Args:
        settings (Kratos.Parameters): The Kratos parameters containing the settings for the Flux wrapper.
        model (Kratos.Model): The Kratos model containing the problem's data.
        solver_name (str): The name of the solver associated with this wrapper.

    Returns:
        FluxWrapper: An instance of the FluxWrapper class.
    """
    return FluxWrapper(settings, model, solver_name)

class FluxWrapper(CoSimulationSolverWrapper):
    """
    This class provides a dedicated Kratos wrapper for Flux 2D and 3D solvers.

    This wrapper inherits from CoSimulationSolverWrapper and defines the standard Kratos wrapper functions for Flux solvers.

    Attributes:
        settings (Kratos.Parameters): The Kratos parameters containing the settings for the Flux wrapper.
        model (Kratos.Model): The Kratos model containing the problem's data.
        solver_name (str): The name of the solver associated with this wrapper.

    Methods:
        __init__(self, settings, model, solver_name): Initialize the Flux wrapper instance.

    Note:
        The FluxWrapper class serves as a bridge between the Flux solver and the Kratos CoSimulation framework.
    """

    def __init__(self, settings, model, solver_name):
        """
        Initialize the FluxWrapper instance.

        Args:
            settings (Kratos.Parameters): The Kratos parameters containing the settings for the Flux wrapper.
            model (Kratos.Model): The Kratos model containing the problem's data.
            solver_name (str): The name of the solver associated with this wrapper.
        """
        super().__init__(settings, model, solver_name)

        # Check for a subprocess call
        world_rank = 0
        if KM.IsDistributedRun() :
            world_data_comm = KM.ParallelEnvironment.GetDataCommunicator("World")
            world_rank = world_data_comm.Rank()

        # Prevent any launch of a FluxServer from a subprocess
        if world_rank != 0:
            return

        # Set default settings and validate JSON settings
        self.settings["solver_wrapper_settings"].ValidateAndAssignDefaults(self.GetDefaultSolverWrapperParameters())

        # ---------------------------
        # Import FluxPythonServer API
        # ---------------------------
        # Getting the path to the FluxPythonServer API
        flux_install_path = self.settings["solver_wrapper_settings"]["flux_install_path"].GetString()

        # If not defined, importing the local version of FluxPythonServer API
        if not flux_install_path:
            current_file_path = Path(__file__).resolve()
            flux_install_path = current_file_path.parent.parent / "Api" / "Python"
        sys.path.append(str(flux_install_path))

        # Importing the FluxPythonServer API
        import FluxPythonServer as FXapi

        # Revert path
        sys.path.remove(str(flux_install_path))

        # ---------------------
        # Start the Flux server
        # ---------------------

        # Get arguments from JSON file
        application = self.settings["solver_wrapper_settings"]["application"].GetString()
        working_directory = self.settings["solver_wrapper_settings"]["working_directory"].GetString()
        # execution_mode = self.settings["solver_wrapper_settings"]["execution_mode"].GetString()
        echo_level = self.settings["solver_wrapper_settings"]["echo_level"].GetInt()
        self.extruded_layers = self.settings["solver_wrapper_settings"]["extruded_layers"].GetInt()
        self.full_device = self.settings["solver_wrapper_settings"]["rebuilt_full_device"].GetString().upper() == "YES" 

        # Start the Flux server and store the instance
        cs_tools.cs_print_info(f"Working directory : {working_directory}")
        cs_tools.cs_print_info(f"Application : {application}")
        cs_tools.cs_print_info(f"Name  : {solver_name}")
        self.flux = FXapi.FluxServer(working_directory, application=application, mode="debug", langage="eng", name=solver_name, verbose=echo_level>0)

        # ---------------------
        # Open the Flux project
        # ---------------------

        # Get argument from JSON file
        flux_project = self.settings["solver_wrapper_settings"]["flux_project"].GetString()
        self.flux.LoadProject(flux_project)

        # Check Json project settings
        self.CheckFluxProject(solver_name)

        # Check physics
        self.flux.CheckPhysics(raise_error=True)

        # -----------------------------
        # Open the MultiPhysics session
        # -----------------------------

        # Get arguments from JSON file
        scenario = self.settings["solver_wrapper_settings"]["scenario"].GetString()
        solved_flux_project = "KratosSolved_"+flux_project

        # Save the project with a new name to prevent any initial project crushing
        self.flux.SaveProjectAs(solved_flux_project)

        # Open a MultiPhysics session
        self.flux.OpenMultiPhysics(scenario)

        # TODO : Should execute a dedicated checker before creating supports, collections and exports

        # ---------------------------------------------------------------------
        # Creation of the model parts and its associates variables in the Kratos environnement
        # --------------------------------------------------------------------------------

        # self.settings["solver_wrapper_settings"].ValidateAndAssignDefaults(settings_defaults)
        model_part_name = self.settings["solver_wrapper_settings"]["main_model_part_name"].GetString()
        model_part_utilities.CreateMainModelPartsFromCouplingDataSettings(self.settings["data"], self.model, self.name)

        # If historical values are required
        model_part_utilities.AllocateHistoricalVariablesFromCouplingDataSettings(self.settings["data"], self.model, self.name)

        # ---------------------------------------------------------------------
        # Creation of the required DataExport instances for each Kratos imports
        # ---------------------------------------------------------------------
        # For each Kratos imported data
        for import_data_name in self.settings["solver_wrapper_settings"]["import_data"].GetStringArray() :
            # Check additional info
            self.settings["data"][import_data_name]["additional_info"].ValidateAndAssignDefaults(self.GetDefaultDataAdditionalInfoParameters())

            # Get arguments of 'data' from JSON file
            model_part_name = self.settings["data"][import_data_name]["model_part_name"].GetString()
            flux_quantity = self.settings["data"][import_data_name]["additional_info"]["flux_quantity"].GetString()
            flux_region_list = self.settings["data"][import_data_name]["additional_info"]["flux_region_list"].GetStringArray()
            initial_value = self.settings["data"][import_data_name]["additional_info"]["initial_value"].GetDouble()
            variable_name = self.settings["data"][import_data_name]["variable_name"].GetString()
            location = self.settings["data"][import_data_name]["location"].GetString()

            # Get general argument
            scenario_type = self.settings["solver_wrapper_settings"]["scenario_type"].GetString()

            # Check entities
            self.CheckDataEntities (import_data_name)

            # Create the Flux DataSupport instance
            data_support = self.flux.CreateDataSupportOnRegion(model_part_name, flux_region_list, self.extruded_layers, self.full_device)

            # Set averaging interval in the case of an instantaneous or a single step scenario type
            if scenario_type == "time_inst" or scenario_type == "single_step":
                # Create the DataExport instance
                self.flux.CreateDataExport (import_data_name,data_support,flux_quantity)

            # Set averaging interval in the case of a time-averaged scenario type
            elif scenario_type == "time_avg" :
                # Get arguments from JSON file
                avg_start_time = self.settings["solver_wrapper_settings"]["avg_start_time"].GetDouble()
                avg_end_time = self.settings["solver_wrapper_settings"]["avg_end_time"].GetDouble()
                interval = [avg_start_time , avg_end_time]
                # Create the DataExport instance
                self.flux.CreateDataExport (import_data_name,data_support,flux_quantity,interval)

            # Import nodes and elements of the model part in the Kratos Environnement
            interface_config = {
                    "model_part_name"  : model_part_name,
                    "initial_value"    : initial_value,
                    "variable_name"    : variable_name,
                    "location"         : location
                    }
            self.ImportCouplingInterface(interface_config)

        # ---------------------------------------------------------------------
        # Creation of the required DataImport instances for each Kratos exports
        # ---------------------------------------------------------------------
        # For each Kratos exported data
        for export_data_name in self.settings["solver_wrapper_settings"]["export_data"].GetStringArray():
            # Check additional info
            self.settings["data"][export_data_name]["additional_info"].ValidateAndAssignDefaults(self.GetDefaultDataAdditionalInfoParameters())

            # Get arguments from JSON file
            model_part_name = self.settings["data"][export_data_name]["model_part_name"].GetString()
            flux_quantity = self.settings["data"][export_data_name]["additional_info"]["flux_quantity"].GetString()
            flux_region_list = self.settings["data"][export_data_name]["additional_info"]["flux_region_list"].GetStringArray()
            initial_value = self.settings["data"][export_data_name]["additional_info"]["initial_value"].GetDouble()
            variable_name = self.settings["data"][export_data_name]["variable_name"].GetString()
            location = self.settings["data"][export_data_name]["location"].GetString()

            # Check entities
            self.CheckDataEntities (export_data_name)

            # Create the Flux DataSupport instance
            data_support = self.flux.CreateDataSupportOnRegion(model_part_name,flux_region_list,self.extruded_layers,self.full_device)

            # Import nodes and elements of the model part in the Kratos Environnement
            interface_config = {
                    "model_part_name"  : model_part_name,
                    "initial_value"    : initial_value,
                    "variable_name"    : variable_name,
                    "location"         : location
                    }
            self.ImportCouplingInterface(interface_config)

            # Create the DataImport instance
            self.flux.CreateDataImport (export_data_name,data_support,flux_quantity)

        # Force first step initialize
        super().Initialize()

    def __del__(self):
        """
        Destructor for the FluxWrapper instance.

        This method is responsible for cleaning up resources associated with the FluxWrapper instance before it is deleted.

        Args:
            self: An instance of the FluxWrapper class.

        Note:
            If the instance has a 'flux' attribute, this method closes the MultiPhysics session, saves and closes the Flux project.
        """

        if hasattr(self, 'flux') :
            # Close the MultiPhysic session
            self.flux.CloseMultiPhysics()
            # Save and close the Flux project
            self.flux.SaveProject()
            self.flux.CloseProject()

    def CheckFluxProject (self,solver_name):
        """
        Check the compatibility of the Flux project with the defined scenario.

        This function checks whether the Flux project defined in the JSON file exists and if the scenario defined in the JSON file is compatible with the Flux project.

        Args:
            self: An instance of the FluxWrapper class.
            solver_name: The name of the solver.

        Note:
            This function verifies that the required Flux project is set up correctly for the given scenario.
        """

        # Check Scenario
        scenario_name = self.settings["solver_wrapper_settings"]["scenario"].GetString()
        if not self.flux.IsExistInstance('Scenario',scenario_name) :
            raise Exception ("Wrong JSON setting : The scenario "+scenario_name+" defined for solver "+solver_name+" does not exist in Flux project.")

        # Check scenario type
        scenario_type = self.settings["solver_wrapper_settings"]["scenario_type"].GetString()
        if scenario_type != "single_step" and scenario_type != "time_inst" and scenario_type != "time_avg" :
            raise Exception ("Wrong JSON setting : The scenario type " + scenario_type + " defined for solver " + solver_name + " must be single_step, time_inst or solver.")

        # Check compatibility with project application
        application = self.flux.GetApplication()
        if application.find("Transient") >= 0 and scenario_type.find("time") <0:
            raise Exception ("Wrong JSON setting : The scenario_type must be time_inst or time_avg when solving  " + application + " Flux application.")
        if application.find("Transient") < 0 and scenario_type.find("single_step") < 0:
            raise Exception ("Wrong JSON setting : The scenario_type must be single_step when solving " + application + " Flux application.")

    def CheckDataEntities(self, data_name) :
        """
        Check the existence of entities defined in the JSON file within the Flux project.

        This function verifies whether the entities specified in the JSON file exist within the Flux project. It checks for the correctness of flux_quantity and the presence of regions in the project.

        Args:
            self: An instance of the FluxWrapper class.
            data_name: The name of the data.

        Raises:
            Exception: If the specified flux_quantity or region does not exist in the Flux project.

        Note:
            This function ensures that the data entities required for the CoSimulation process are properly configured in the Flux project.
        """

        # Check flux_quantity
        flux_quantity = self.settings["data"][data_name]["additional_info"]["flux_quantity"].GetString()
        if not flux_quantity in ["JOULE_LOSSES","IRON_LOSSES","TKELVIN","TCELSIUS"] :
            if not self.flux.IsExistInstance('SpatialParameter',flux_quantity) :
                raise Exception ("Wrong JSON setting : The spatial parameter "+flux_quantity+" defined in data "+data_name+" does not exist in Flux project.")

        # Check regions
        dimension = self.flux.GetDimension()
        if dimension == 2 : region_type = "RegionFace"
        if dimension == 3 : region_type = "RegionVolume"
        for region in self.settings["data"][data_name]["additional_info"]["flux_region_list"].GetStringArray() :
            if not self.flux.IsExistInstance(region_type,region):
                raise Exception ("Wrong JSON setting : The region " + region + " defined in data " + data_name + " does not exist in Flux project.")


    def Initialize(self):
        """
        Initialize the CoSimulation step.

        This method initializes the CoSimulation step by calling the corresponding method in the Flux solver.

        Args:
            self: An instance of the FluxWrapper class.
        """

        self.flux.InitializeStep()

    def Finalize(self):
        """
        Finalize the CoSimulation step.

        This method finalizes the CoSimulation step by calling the corresponding method in the Flux solver.

        Args:
            self: An instance of the FluxWrapper class.
        """

        self.flux.FinalizeStep()

    def AdvanceInTime(self, current_time):
        """
        Advance the solution in time based on the specified scenario type.

        This function advances the solution in time based on the scenario type defined in the JSON settings.

        Args:
            self: An instance of the FluxWrapper class.
            current_time: The current time.

        Returns:
            float: The new time after the time advancement.

        Note:
            The behavior of this method depends on the scenario type specified in the settings.
        """

        # In  the case of current_time = -1.0
        if current_time < 0.0 :
            return 0.0

        # Get the Flux scenario type from JSON file
        scenario_type = self.settings["solver_wrapper_settings"]["scenario_type"].GetString()
        if scenario_type == "single_step" :
            # Nothing to do in the case
            return 0.0
        elif scenario_type == "time_avg" :
            # Delete all solved steps to reactivate first step of the scenario
            self.flux.DeleteAllSolvedSteps()
            # No pertinent active time
            return 0.0
        elif scenario_type == "time_inst" :
            self.flux.ActivateNextStep()
            self.flux.InitializeStep()
            new_time = self.flux.GetCurrentTime()
            if current_time == None : return 0.0
            else : return new_time

    def SolveSolutionStep(self):
        """
        Solves the current time step including export/import processes.

        Args:
            self: An instance of the class.
        """

        # ------------------------
        # Solve the required steps
        # ------------------------
        # Get scenario type from JSON file
        scenario_type = self.settings["solver_wrapper_settings"]["scenario_type"].GetString()

        # In the case of a single step or an instantaneous scenario
        if scenario_type == "single_step" or scenario_type == "time_inst" :
            # Export all data from Kratos to Flux
            self.ExportAllData()
            # Solve the current step
            self.flux.SolveStep()
            # Import all data from Flux to Kratos
            self.ImportAllData()

        # In the case of a time averaged scenario
        elif scenario_type  == "time_avg" :
            # If scenario is already solved
            if not self.flux.HasNextStep() :
                self.flux.DeleteAllSolvedSteps()
            # Export all data from Kratos to Flux
            self.ExportAllData()
            # Solve all steps
            self.flux.SolveAllSteps()
            # Import all data from Flux to Kratos
            self.ImportAllData()

    def _GetIOType(self):
        """
        Get the Flux I/O type dedicated to handling import/export data.

        Returns:
            str: The string "flux_io".

        Note:
            This method specifies the dedicated I/O type for handling data import and export using Flux.

        """

        return "flux_io"

    def _GetDataCommunicator(self):
        """
        Get the data communicator dedicated to handling import/export data.

        Returns:
            DataCommunicator: The data communicator instance.

        Note:
            This method returns the data communicator to be used for data import and export.
            In this case, the solver does not support MPI, so the RankZeroDataCommunicator is returned.

        """

        return data_communicator_utilities.GetRankZeroDataCommunicator()

    def __RunExecutable(self):
        """
        Run the executable associated with the solver.

        Note:
            This function is not implemented and raises a NotImplementedError.
            The Flux dedicated wrapper does not currently support running an executable.

        Args:
            self: An instance of the FluxWrapper class.

        Raises:
            NotImplementedError: If the method is called, indicating that running an executable is not supported.

        """

        # Not implemented
        raise NotImplementedError('RunExecutable is not implemented for the Flux dedicated wrapper')

    def ImportCouplingInterface(self, interface_config): 
        """
        Create a Kratos model-part from a Flux DataSupport.

        This function creates a Kratos model-part by importing data from a Flux DataSupport based on the provided interface configuration.

        Args:
            self: An instance of the FluxWrapper class.
            interface_config (dict): A dictionary containing the interface configuration.

        Note:
            The method extracts information from the interface_config dictionary, retrieves data from Flux DataSupport, and creates Kratos nodes and elements in the model part.

        """

        # Get arguments
        model_part_name = interface_config["model_part_name"]
        initial_value = interface_config["initial_value"]
        variable_name = interface_config["variable_name"]
        location = interface_config["location"]

        # Get the DataSupport instance by its name
        data_support = self.flux.GetDataSupportByName(model_part_name)

        # Get the mesh of the Data support
        [node_id,node_coordinates,element_id,element_type,element_part,element_nodes] = data_support.GetMesh(connectivity='KRATOS')

        # Get the model part by its name
        model_part = self.model.GetModelPart(model_part_name)
        model_part.AddProperties(KM.Properties(1))

        # Get the Kratos variable by its name
        variable = KM.KratosGlobals.GetVariable(variable_name)

        # TODO: Clean all existing nodes and elements of the model part ?
        # Create the nodes and affect initial value
        for i in range(len(node_id)) :
            node = model_part.CreateNewNode(node_id[i],node_coordinates[i][0],node_coordinates[i][1],node_coordinates[i][2])
            if location == 'node_historical' :
                node.SetSolutionStepValue(variable, 0, initial_value)
            else :
                node.SetValue(variable, initial_value)
        # Create the elements
        for i in range(len(element_id)) :
            model_part.CreateNewElement(element_type[i],element_id[i],element_nodes[i],model_part.GetProperties()[1])

    def ExportCouplingInterface(self, interface_config):
        """
        Placeholder method for transferring a mesh from Kratos to Flux.

        Note:
            This function is not needed at the moment and is provided as a placeholder.

        Args:
            interface_config: Configuration data for the coupling interface.

        """

        pass

    def ImportData(self, data_config):
        """
        Transfer specified data from Flux to Kratos.

        This function transfers specified data from the Flux solver to the Kratos model based on the provided data configuration.

        Args:
            self: An instance of the FluxWrapper class.
            data_config (dict): A dictionary containing the data configuration.

        Note:
            The method extracts information from the data_config dictionary, retrieves data from Flux, and updates the Kratos model part nodes.

        """

        # Get arguments
        import_data_name = data_config["import_data_name"]
        model_part_name = data_config["model_part_name"]
        variable_name = data_config["variable_name"]
        location =  data_config["location"]

        # Get the DataExport instance by its name
        data_export = self.flux.GetDataExportByName(import_data_name)

        # Get the model part by its name
        model_part = self.model.GetModelPart(model_part_name)

        # Get the Kratos variable by its name
        variable = KM.KratosGlobals.GetVariable(variable_name)

        # Get the values at nodes
        [node_ids,node_values] = data_export.GetData()

        # Affect values to nodes of the Kratos model part
        for node in model_part.Nodes :
            if location == 'node_historical' :
                node.SetSolutionStepValue(variable, 0, node_values[node_ids.index(node.Id)])
            else :
                node.SetValue(variable, node_values[node_ids.index(node.Id)])

    def ImportAllData(self) :
        """
        Import all specified data from Flux into Kratos.

        This function imports all the specified data from the Flux solver into the Kratos model based on the defined data configurations.

        Args:
            self: An instance of the FluxWrapper class.

        Note:
            The method iterates over the list of importable data defined in the JSON settings and uses the ImportData method for each data entry.

        """

        for import_data in self.settings["solver_wrapper_settings"]["import_data"].GetStringArray():
            data_config = {
                "import_data_name" : import_data,
                "model_part_name"  : self.settings["data"][import_data]["model_part_name"].GetString(),
                "variable_name"    : self.settings["data"][import_data]["variable_name"].GetString(),
                "location"         : self.settings["data"][import_data]["location"].GetString()
            }
            self.ImportData(data_config)

    def ExportData(self, data_config):
        """
        Transfer specified data from Kratos to Flux.

        This function transfers specified data from the Kratos model to the Flux solver based on the provided data configuration.

        Args:
            self: An instance of the FluxWrapper class.
            data_config (dict): A dictionary containing the data configuration.

        Note:
            The method extracts information from the data_config dictionary, retrieves data from Kratos, and updates the Flux DataImport instance.

        """

        # Check if method is called from
        # GaussSeidelStrongCoupledSolver.__CommunicateIfTimeStepNeedsToBeRepeated()
        if "type" in data_config :
            # In the case, what should I do ?
            # If nothing, it seems to work
            return

        # Get arguments
        export_data_name = data_config["export_data_name"]
        model_part_name = data_config["model_part_name"]
        variable_name = data_config["variable_name"]
        location =  data_config["location"]

        # Get the DataImport instance by its name
        data_import = self.flux.GetDataImportByName(export_data_name)

        # Get the model part by its name
        model_part = self.model.GetModelPart(model_part_name)

        # Get the Kratos variable by its name
        variable = KM.KratosGlobals.GetVariable(variable_name)

        # Get values at nodes from the Kratos model part
        node_ids = []
        node_values = []
        for node in model_part.Nodes:
            node_ids.append(node.Id)
            if location == 'node_historical' :
                node_values.append(node.GetSolutionStepValue(variable, 0))
            else :
                node_values.append(node.GetValue(variable))

        # Update current values of the DataImport instance
        data_import.SetData(node_ids,node_values)

    def ExportAllData (self) :
        """
        Export all specified data from Kratos to Flux.

        This function exports all the specified data from the Kratos model to the Flux solver. The data to be exported is defined in the JSON settings.

        Args:
            self: An instance of the FluxWrapper class.

        Note:
            The method iterates over the list of exportable data defined in the JSON settings and uses the ExportData method for each data entry.

        """

        # Export all data from Kratos to Flux and update data in Flux
        for export_data in self.settings["solver_wrapper_settings"]["export_data"].GetStringArray():
            data_config = {
                "export_data_name" : export_data,
                "model_part_name"  : self.settings["data"][export_data]["model_part_name"].GetString(),
                "variable_name"    : self.settings["data"][export_data]["variable_name"].GetString(),
                "location"         : self.settings["data"][export_data]["location"].GetString()
            }
            self.ExportData(data_config)

    @staticmethod
    def GetDefaultSolverWrapperParameters():
        """
        Get the default parameters for configuring the solver wrapper.

        Returns:
            Kratos.Parameters: The default parameters for configuring the solver wrapper.

        Note:
            These default parameters define the basic settings needed for setting up the Flux solver wrapper.
        """

        return KM.Parameters("""{
            "main_model_part_name" : "",
            "application"          : "Flux3D",
            "extruded_layers"      : 0,
            "rebuilt_full_device"  : "yes",
            "flux_project"         : "",
            "execution_mode"       : "RELEASE",
            "working_directory"    : ".",
            "scenario"             : "",
            "scenario_type"        : "single_step",
            "avg_start_time"       : 0.0,
            "avg_end_time"         : 0.0,
            "echo_level"           : 0,
            "export_data"          : [ ],
            "import_data"          : [ ],
            "flux_install_path"    : ""
        }""")

    @staticmethod
    def GetDefaultDataAdditionalInfoParameters():
        """
        Get the default parameters for data additional information.

        Returns:
            Kratos.Parameters: The default parameters for data additional information.

        Note:
            These default parameters define the additional information associated with data settings in the Flux solver.
        """

        return KM.Parameters("""{
            "initial_value"     : 0.0,
            "flux_quantity"     : "",
            "flux_region_list"  : [ ]
        }""")