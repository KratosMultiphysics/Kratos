# Importing the Kratos Library
import KratosMultiphysics as KM
import subprocess

# Importing the base class
from KratosMultiphysics.CoSimulationApplication.base_classes.co_simulation_solver_wrapper import CoSimulationSolverWrapper

# Other imports
from KratosMultiphysics.CoSimulationApplication.utilities import model_part_utilities
import KratosMultiphysics.CoSimulationApplication.co_simulation_tools as cs_tools
from KratosMultiphysics.vtk_output_process import VtkOutputProcess
import time
import sys
import os
from sys import platform

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
        
        self.name = solver_name
        self.PrintLog("Initiliazing AcuSolve Wrapper...")
        super().__init__(settings, model, solver_name)
        # Set default settings and validate JSON settings
        solver_wrapper_settings_defaults = KM.Parameters("""{
            "main_model_part_name" : "",
            "application" : "AcuSolve",
            "problem" : "",
            "input_file" : "",
            "working_directory" : "",
            "problem_directory" : ".",
            "num_processors" : 1,
            "num_threads" : 1,
            "time_increment" : 0.0,
            "region" : "",
            "gpu_flag" : "FALSE",
            "restart_flag" : "FALSE",
            "fast_restart_flag" : "FALSE",
            "echo_level" : 1,
            "export_data"      : [ ],
            "import_data"      : [ ],
            "import_meshes"      : [ ],
            "execution_mode"  : "RELEASE"
        }""")

        self.settings["solver_wrapper_settings"].ValidateAndAssignDefaults(solver_wrapper_settings_defaults)
        
        # --------------------------------------------------------------------------------
        # Creation of the model parts and its associates variables in the Kratos environnement
        # --------------------------------------------------------------------------------
        model_part_utilities.CreateMainModelPartsFromCouplingDataSettings(self.settings["data"], self.model, self.name)
        model_part_utilities.AllocateHistoricalVariablesFromCouplingDataSettings(self.settings["data"], self.model, self.name)        

        # ---------------------
        # Get arguments from JSON file
        application         = self.settings["solver_wrapper_settings"]["application"].GetString()
        working_directory   = self.settings["solver_wrapper_settings"]["working_directory"].GetString()
        problem_directory   = self.settings["solver_wrapper_settings"]["problem_directory"].GetString()
        region              = self.settings["solver_wrapper_settings"]["region"].GetString()
        echo_level          = self.settings["solver_wrapper_settings"]["echo_level"].GetInt()
        problem             = self.settings["solver_wrapper_settings"]["problem"].GetString()
        inputFile           = self.settings["solver_wrapper_settings"]["input_file"].GetString()
        np                  = self.settings["solver_wrapper_settings"]["num_processors"].GetInt()
        nt                  = self.settings["solver_wrapper_settings"]["num_threads"].GetInt()
        gpu                 = self.settings["solver_wrapper_settings"]["gpu_flag"].GetString()
        rst                 = self.settings["solver_wrapper_settings"]["restart_flag"].GetString()
        frst                = self.settings["solver_wrapper_settings"]["restart_flag"].GetString()
        self.deltaT         = self.settings["solver_wrapper_settings"]["time_increment"].GetDouble()
        # Evaluate debug mode indicator
        self.debug_mode = self.settings["solver_wrapper_settings"]["execution_mode"].GetString().upper() == "DEBUG"

        # Set internal process variables
        self.advance_in_time = False

        # To be used for distributed runs adaptation :
        #if KM.IsDistributedRun() :
        #    world_data_comm = KM.ParallelEnvironment.GetDataCommunicator("World")
        #    world_rank = world_data_comm.Rank()
        #    if ( world_rank == 0 ) :

        # Start AcuSolve process
        self.PrintLog("Starting AcuSolve process...")
        if platform == "linux" or platform == "linux2":
            if (gpu!="FALSE"):
                cmd = "acuRun -inp "+inputFile+" -pb "+problem+" -dir \""+working_directory+"\" -np "+str(np)+" -nt "+str(nt)+" -gpu "+gpu+" -verbose "+str(echo_level)+" &"         
            elif (rst!="FALSE"):
                cmd = "acuRun -inp "+inputFile+" -pb "+problem+" -dir \""+working_directory+" -np "+str(np)+" -nt "+str(nt)+" -rst "+" -verbose "+str(echo_level)+" &"         
            elif (frst!="FALSE"):
                cmd = "acuRun -inp "+inputFile+" -pb "+problem+" -dir \""+working_directory+"\" -np "+str(np)+" -nt "+str(nt)+" -frst "+" -verbose "+str(echo_level)+" &"         
            else:
                cmd = "acuRun -inp "+inputFile+" -pb "+problem+" -dir \""+working_directory+"\" -np "+str(np)+" -nt "+str(nt)+" -verbose "+str(echo_level)+" &"       
            self.PrintLog(f"AcuSolve command : {cmd}")
            os.system(cmd)
        elif platform == "win32":
            if (gpu!="FALSE"):
                cmd = "acuRun.bat -inp "+inputFile+" -pb "+problem+" -dir \""+working_directory+"\" -pdir \""+problem_directory+"\" -np "+str(np)+" -nt "+str(nt)+" -gpu "+gpu+" -verbose "+str(echo_level)+" -lbuff"         
            elif (rst!="FALSE"):
                cmd = "acuRun.bat -inp "+inputFile+" -pb "+problem+" -dir \""+working_directory+"\" -pdir \""+problem_directory+"\" -np "+str(np)+" -nt "+str(nt)+" -rst "+" -verbose "+str(echo_level)+" -lbuff"         
            elif (frst!="FALSE"):
                cmd = "acuRun.bat -inp "+inputFile+" -pb "+problem+" -dir \""+working_directory+"\" -pdir \""+problem_directory+"\" -np "+str(np)+" -nt "+str(nt)+" -frst "+" -verbose "+str(echo_level)+" -lbuff"        
            else:
                cmd = "acuRun.bat -inp "+inputFile+" -pb "+problem+" -dir \""+working_directory+"\" -pdir \""+problem_directory+"\" -np "+str(np)+" -nt "+str(nt)+" -verbose "+str(echo_level)+" -lbuff"  
            self.PrintLog(f"AcuSolve command : {cmd}")
            subprocess.Popen(cmd)

    def Initialize(self):
        # The initialize process must be called first after starting AcuSolve
        # The actions to follow are dicted by the Acusolve process :
        #   1. Establish the connection by using CosimIO
        #   2. Import all meshes from Acusolve
        #   3. Import initiale values of temperature from Acusolve
        #   4. Export initiale values of heat sources to Acusolve
        self.PrintLog("Initialize")

        # Connect Kratos to Acusolve with CosimIO
        super().Initialize()

        # Import all meshes from Acusolve
        for model_part_name in self.settings["solver_wrapper_settings"]["import_meshes"].GetStringArray():
             interface_config = { "model_part_name" : model_part_name }
             self.ImportCouplingInterface(interface_config)

        # In debub mode, start the VTK output process
        if self.debug_mode:
            self.PrintLog("Starting VTK ouput process")
            vtk_output_configuration = KM.Parameters("""{
                    "model_part_name"        : \""""+model_part_name+"""\",
                    "output_sub_model_parts" : false,
                    "nodal_solution_step_data_variables" : ["HEAT_FLUX"]
                }""")
            self.vtk_output = VtkOutputProcess(self.model, vtk_output_configuration)
            self.vtk_output.ExecuteInitialize()
            self.vtk_output.ExecuteBeforeSolutionLoop()
            self.vtk_output.PrintOutput()

        # Import initial temperature values from Acusolve     
        for data_name in self.settings["solver_wrapper_settings"]["import_data"].GetStringArray():
            data_config = { 
                "type" : "coupling_interface_data",
                "interface_data" : self.GetInterfaceData(data_name)
            }
            self.ImportData(data_config)

        # Export initial heat sources to Acusolve
        for data_name in self.settings["solver_wrapper_settings"]["export_data"].GetStringArray():
            data_config = { 
                "type" : "coupling_interface_data",
                "interface_data" : self.GetInterfaceData(data_name)
            }
            self.ExportData(data_config)
       
        # At this point, Acusolve is waiting to solve a step...
        return

    def AdvanceInTime(self, current_time):
        # The AdvanceInTime process does not affect directly Acusolve
        # but the information to advance in time (or not) must be
        # store to be used when SolveSolutionStep is called
        self.PrintLog("AdvanceInTime")
        # Because the coupled_solver instance calls AdvanceInTime first before solving the first step
        # this first call must be skipped for Acusolve which has already the proper
        # step loaded
        if (current_time > 0.0) : 
            self.advance_in_time = True
        return (current_time+self.deltaT)

    def SolveSolutionStep(self):
        # The SolveSolutionStep process process in Acusolve follows a specific order which
        # must be respected. The order is the following:
        #   1. (If required) Activate the next step in Acusolve
        #   2. Solve the step
        #   3. Export the temperature value to Kratos
        #   4. Wait for and Import the heat sources from Kratos
        # The process starts when the AdvanceInTime or RepeatTimeStep signal is sent to Acusolve
        start_time = time.time()
        self.PrintLog("SolveSolutionStep starts...")

        # In debug mode
        if self.debug_mode :
            self.PrintLog("Exporting current I/O data to VTK output")
            self.vtk_output.PrintOutput()
            self.vtk_output.ExecuteFinalizeSolutionStep()

        # Start the solving process in Acusolve
        if self.advance_in_time : 
            self.PrintLog("Activate next step")
            self.__SendControlSignal("AdvanceInTime")
            self.advance_in_time = False
        else:
            self.PrintLog("Repeat step")
            self.__SendControlSignal("RepeatTimeStep")

        # Wait for and import the temperature values from Acusolve
        self.PrintLog("Import temperature from Acusolve")
        for data_name in self.settings["solver_wrapper_settings"]["import_data"].GetStringArray():
            data_config = { 
                "type" : "coupling_interface_data",
                "interface_data" : self.GetInterfaceData(data_name)
            }
            self.ImportData(data_config)

        # Export heat sources to Acusolve (who is waiting for them)
        self.PrintLog("Export heat sources to Acusolve")
        for data_name in self.settings["solver_wrapper_settings"]["export_data"].GetStringArray():
            data_config = { 
                "type" : "coupling_interface_data",
                "interface_data" : self.GetInterfaceData(data_name)
            }
            self.ExportData(data_config)

        self.PrintLog(f"SolveSolutionStep executed. :: {time.time() - start_time} s.")

        return

    def Finalize(self):
        # Inform Acusolve that the overall process is finished.
        # The order will activated finalizing process in Acusolve.
        self.PrintLog("Finalize")
        self.__SendControlSignal("exit")
        # If the VTK output process is running, finalize it
        if self.debug_mode:
            self.PrintLog("Finalizing VTK output process")
            self.vtk_output.ExecuteFinalize()
        super().Finalize()

    def _GetIOType(self):
        return self.settings["io_settings"]["type"].GetString()

    def __SendControlSignal(self, signal, settings=None):
        # The control signal process is used to 
        # inform Acusolve what to do next, it must be :
        #   - RepeatTimeStep
        #   - AdvanceInTime
        #   - exit
        self.PrintLog(f"Sending signal : {signal}")
        if signal not in ["RepeatTimeStep","AdvanceInTime","exit"]:
            raise Exception("Signal {} is not authorized".format(signal))
        data_config = {
            "type"           : "control_signal",
            "control_signal" : signal,
            "settings"       : settings
        }
        # The signal is sent to Acusolve through the CosimIO export function
        self.ExportData(data_config)

    def PrintLog(self,text):
        return cs_tools.cs_print_info(f"{self.name} :",text)
