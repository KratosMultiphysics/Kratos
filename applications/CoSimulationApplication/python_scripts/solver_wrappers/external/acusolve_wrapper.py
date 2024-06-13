# Importing the Kratos Library
import KratosMultiphysics as KM

# Importing the base class
from KratosMultiphysics.CoSimulationApplication.base_classes.co_simulation_solver_wrapper import CoSimulationSolverWrapper
from KratosMultiphysics.vtk_output_process import VtkOutputProcess

# CoSimulation imports
from KratosMultiphysics.CoSimulationApplication.utilities import model_part_utilities
import KratosMultiphysics.CoSimulationApplication.co_simulation_tools as cs_tools

# System imports
import subprocess
import platform

def Create(settings, model, solver_name):
    return acuSolveWrapper(settings, model, solver_name)

class acuSolveWrapper(CoSimulationSolverWrapper):
    """This class serves as a dedicated Kratos wrapper for acuSolve
    """
    # -----------------------------------------------
    # Definition of standard Kratos wrapper functions
    # -----------------------------------------------
    # Create the AcuSolve Wrapper instance
    def __init__(self, settings, model, solver_name):
        """Constructor of the AcuSolve Wrapper
        self        : an instance of the class
        settings    : a Kratos Parameters object containing the settings
        model       : a Kratos Model object containing the model
        solver_name : the name of the solver
        """
        
        super().__init__(settings, model, solver_name)
        # Set default settings and validate JSON settings
        solver_wrapper_settings_defaults = KM.Parameters("""{
            "main_model_part_name"    : "",
            "application"             : "AcuSolve",
            "problem"                 : "",
            "input_file"              : "",
            "working_directory"       : "",
            "problem_directory"       : ".",
            "num_processors"          : 1,
            "num_threads"             : 1,
            "time_increment"          : 0.0,
            "region"                  : "",
            "gpu_flag"                : "FALSE",
            "restart_flag"            : "FALSE",
            "fast_restart_flag"       : "FALSE",
            "echo_level"              : 1,
            "export_data"             : [ ],
            "import_data"             : [ ],
            "import_meshes"           : [ ],
            "post_process_first_step" :  true
        }""")

        # Validate settings
        self.settings["solver_wrapper_settings"].ValidateAndAssignDefaults(solver_wrapper_settings_defaults)

        # Get solver name
        self.name = solver_name
        
        # --------------------------------------------------------------------------------
        # Creation of the model parts and its associates variables in the Kratos environnement
        # --------------------------------------------------------------------------------
        model_part_utilities.CreateMainModelPartsFromCouplingDataSettings(self.settings["data"], self.model, self.name)
        model_part_utilities.AllocateHistoricalVariablesFromCouplingDataSettings(self.settings["data"], self.model, self.name)        
        cs_tools.cs_print_info(self.name + ": " +  "Run AcuSolve")

        # ---------------------
        # Get arguments from JSON file
        #application                  = self.settings["solver_wrapper_settings"]["application"].GetString()
        working_directory            = self.settings["solver_wrapper_settings"]["working_directory"].GetString()
        problem_directory            = self.settings["solver_wrapper_settings"]["problem_directory"].GetString()
        #region                       = self.settings["solver_wrapper_settings"]["region"].GetString()
        echo_level                   = self.settings["solver_wrapper_settings"]["echo_level"].GetInt()
        problem                      = self.settings["solver_wrapper_settings"]["problem"].GetString()
        inputFile                    = self.settings["solver_wrapper_settings"]["input_file"].GetString()
        np                           = self.settings["solver_wrapper_settings"]["num_processors"].GetInt()
        nt                           = self.settings["solver_wrapper_settings"]["num_threads"].GetInt()
        gpu                          = self.settings["solver_wrapper_settings"]["gpu_flag"].GetString()
        rst                          = self.settings["solver_wrapper_settings"]["restart_flag"].GetString()
        frst                         = self.settings["solver_wrapper_settings"]["restart_flag"].GetString()
        self.deltaT                  = self.settings["solver_wrapper_settings"]["time_increment"].GetDouble()
        self.post_process_first_step = self.settings["solver_wrapper_settings"]["post_process_first_step"].GetBool()

        # ---------------------
        # Launch AcuSolve
        # ---------------------
        common_cmd = inputFile + " -pb " + problem +" -dir "+ working_directory + " -pdir " + problem_directory + " -np "+ str(np) + " -nt " + str(nt)
        verbose_cmd = " -verbose " + str(echo_level)
        platform_details = {
            "Linux"   : ("acuRun -inp ", " &"),
            "Windows" : ("acuRun.bat -inp ", " -lbuff")
        }
        current_platform = platform.system()
        if current_platform in platform_details:
            base_cmd, buffer_cmd = platform_details[current_platform]
            
            if gpu != "FALSE":
                cmd = f"{base_cmd}{common_cmd} -gpu {gpu}{verbose_cmd}{buffer_cmd}"
            elif rst != "FALSE":
                cmd = f"{base_cmd}{common_cmd} -rst {verbose_cmd}{buffer_cmd}"
            elif frst != "FALSE":
                cmd = f"{base_cmd}{common_cmd} -frst {verbose_cmd}{buffer_cmd}"
            else:
                cmd = f"{base_cmd}{common_cmd}{verbose_cmd}{buffer_cmd}"
        else:
            raise Exception("Unsupported operating system detected.")
        cs_tools.cs_print_info(self.name + ": " +  cmd)
        subprocess.run(cmd, shell=True)

    def Initialize(self):
        """ This function initializes the AFS Wrapper
        self : an instance of the class
        """
        super().Initialize()
        for model_part_name in self.settings["solver_wrapper_settings"]["import_meshes"].GetStringArray():
             interface_config = { "model_part_name" : model_part_name }
             self.ImportCouplingInterface(interface_config)
        
        # Post-process if required
        if self.post_process_first_step:
            cs_tools.cs_print_info(self.name + ": " +  "WRITING VTK OUTPUT...........")
            vtk_output_configuration = KM.Parameters("""{
                    "model_part_name"        : \""""+model_part_name+"""\",
                    "output_sub_model_parts" : false,
                    "nodal_solution_step_data_variables" : ["HEAT_FLUX"]
                }""")
            self.vtk_output = VtkOutputProcess(self.model, vtk_output_configuration)
            self.vtk_output.ExecuteInitialize()
            self.vtk_output.ExecuteBeforeSolutionLoop()
            self.vtk_output.PrintOutput()

    def ExportData(self, data_config):
        """This function exports data to the AFS Wrapper
        self        : an instance of the class
        data_config : a dictionary containing the data configuration
        """
        if data_config["type"] == "repeat_time_step" and data_config["repeat_time_step"] == True:
            cs_tools.cs_print_info(self.name + ": " +  "control signal : RepeatTimeStep")
            self.__SendControlSignal("RepeatTimeStep")
            #self.__SendControlSignal("Repeat")
            return # we control the ext solver, no need for sending a repeat_time_step signal
        elif data_config["type"] == "repeat_time_step" and data_config["repeat_time_step"] == False:
           return
        #    cs_tools.cs_print_info(self.name + ": " +  "control signal : NoRepeatTimeStep")
        #    self.__SendControlSignal("NoRepeatTimeStep")  
        super().ExportData(data_config)

    def AdvanceInTime(self, current_time):
        """This function advances the solution in time
        self        : an instance of the class
        current_time: the current time
        """
        if current_time > 0.0:
            cs_tools.cs_print_info(self.name + ": " +  "control signal : AdvanceInTime")
            self.__SendControlSignal("AdvanceInTime")
            #self.__SendControlSignal("Advance")
        cs_tools.cs_print_info(self.name + ": " +  self.deltaT)
        return (current_time+self.deltaT)  # TODO find a better solution here... maybe get time from solver through IO
        #return 0.0

    def SolveSolutionStep(self):
        """This function solves the solution step
        self : an instance of the class
        """
        for data_name in self.settings["solver_wrapper_settings"]["export_data"].GetStringArray():
            data_config = { 
                "type" : "coupling_interface_data",
                "interface_data" : self.GetInterfaceData(data_name)
            }
            self.ExportData(data_config)

        super().SolveSolutionStep()

        for data_name in self.settings["solver_wrapper_settings"]["import_data"].GetStringArray():
            data_config = { 
                "type" : "coupling_interface_data",
                "interface_data" : self.GetInterfaceData(data_name)
            }
            self.ImportData(data_config)

    def Finalize(self):
        """Finalization of the coupled solver.
        self : The coupled solver object
        """
        cs_tools.cs_print_info(self.name + ": " +  "control signal : exit")
        self.__SendControlSignal("exit")
        self.vtk_output.ExecuteFinalize()
        super().Finalize()

    def _GetIOType(self):
        return self.settings["io_settings"]["type"].GetString()

    def __SendControlSignal(self, signal, settings=None):
        data_config = {
            "type"           : "control_signal",
            "control_signal" : signal,
            "settings"       : settings
        }
        self.ExportData(data_config)