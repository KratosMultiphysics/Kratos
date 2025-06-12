# Importing the Kratos Library
import KratosMultiphysics as KM
import subprocess

# Importing the base class
from KratosMultiphysics.CoSimulationApplication.base_classes.co_simulation_solver_wrapper import CoSimulationSolverWrapper

# Other imports
from KratosMultiphysics.CoSimulationApplication.utilities import model_part_utilities
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
            "import_meshes"      : [ ]
        
        }""")

        self.settings["solver_wrapper_settings"].ValidateAndAssignDefaults(solver_wrapper_settings_defaults)
        
        # --------------------------------------------------------------------------------
        # Creation of the model parts and its associates variables in the Kratos environnement
        # --------------------------------------------------------------------------------
        model_part_utilities.CreateMainModelPartsFromCouplingDataSettings(self.settings["data"], self.model, self.name)
        model_part_utilities.AllocateHistoricalVariablesFromCouplingDataSettings(self.settings["data"], self.model, self.name)        
        print ("Run AcuSolve")
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

        #if KM.IsDistributedRun() :
        #    world_data_comm = KM.ParallelEnvironment.GetDataCommunicator("World")
        #    world_rank = world_data_comm.Rank()

        #    if ( world_rank == 0 ) :

                # ---------------------
                # lauch AcuSolve
                # ---------------------

        if platform == "linux" or platform == "linux2":
            if (gpu!="FALSE"):
                cmd = "acuRun -inp "+inputFile+" -pb "+problem+" -dir "+working_directory+" -pdir "+problem_directory+" -np "+str(np)+" -nt "+str(nt)+" -gpu "+gpu+" -verbose "+str(echo_level)+" &"         
            elif (rst!="FALSE"):
                cmd = "acuRun -inp "+inputFile+" -pb "+problem+" -dir "+working_directory+" -pdir "+problem_directory+" -np "+str(np)+" -nt "+str(nt)+" -rst "+" -verbose "+str(echo_level)+" &"         
            elif (frst!="FALSE"):
                cmd = "acuRun -inp "+inputFile+" -pb "+problem+" -dir "+working_directory+" -pdir "+problem_directory+" -np "+str(np)+" -nt "+str(nt)+" -frst "+" -verbose "+str(echo_level)+" &"         
            else:
                cmd = "acuRun -inp "+inputFile+" -pb "+problem+" -dir "+working_directory+" -pdir "+problem_directory+" -np "+str(np)+" -nt "+str(nt)+" -verbose "+str(echo_level)+" &"       
            print (cmd)
            os.system(cmd)
        elif platform == "win32":
            if (gpu!="FALSE"):
                cmd = "acuRun.bat -inp "+inputFile+" -pb "+problem+" -dir "+working_directory+" -pdir "+problem_directory+" -np "+str(np)+" -nt "+str(nt)+" -gpu "+gpu+" -verbose "+str(echo_level)+" -lbuff"         
            elif (rst!="FALSE"):
                cmd = "acuRun.bat -inp "+inputFile+" -pb "+problem+" -dir "+working_directory+" -pdir "+problem_directory+" -np "+str(np)+" -nt "+str(nt)+" -rst "+" -verbose "+str(echo_level)+" -lbuff"         
            elif (frst!="FALSE"):
                cmd = "acuRun.bat -inp "+inputFile+" -pb "+problem+" -dir "+working_directory+" -pdir "+problem_directory+" -np "+str(np)+" -nt "+str(nt)+" -frst "+" -verbose "+str(echo_level)+" -lbuff"        
            else:
                cmd = "acuRun.bat -inp "+inputFile+" -pb "+problem+" -dir "+working_directory+" -pdir "+problem_directory+" -np "+str(np)+" -nt "+str(nt)+" -verbose "+str(echo_level)+" -lbuff"  
            print (cmd)
            subprocess.Popen(cmd)


    def Initialize(self):
        super().Initialize()
        for model_part_name in self.settings["solver_wrapper_settings"]["import_meshes"].GetStringArray():
             interface_config = { "model_part_name" : model_part_name }
             self.ImportCouplingInterface(interface_config)
        
        print("WRITING VTK OUTPUT...........")
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
        if data_config["type"] == "repeat_time_step" and data_config["repeat_time_step"] == True:
            print("control signal : RepeatTimeStep")
            self.__SendControlSignal("RepeatTimeStep")
            #self.__SendControlSignal("Repeat")
            return # we control the ext solver, no need for sending a repeat_time_step signal
        elif data_config["type"] == "repeat_time_step" and data_config["repeat_time_step"] == False:
           return
        #    print("control signal : NoRepeatTimeStep")
        #    self.__SendControlSignal("NoRepeatTimeStep")  
        super().ExportData(data_config)

    def AdvanceInTime(self, current_time):
        if (current_time > 0.0):
            print("control signal : AdvanceInTime")
            self.__SendControlSignal("AdvanceInTime")
            #self.__SendControlSignal("Advance")
        print(self.deltaT) 
        return (current_time+self.deltaT)  # TODO find a better solution here... maybe get time from solver through IO
        #return 0.0

    def SolveSolutionStep(self):
        #for model_part_name in self.settings["solver_wrapper_settings"]["import_meshes"].GetStringArray():
            #model_part = self.model.GetModelPart(model_part_name)
            #print (mode_part.Name())
            #print("PRINTING VTK OUTPUT...................................",self.name,flush = True)
            #self.vtk_output.PrintOutput()
            #self.vtk_output.ExecuteFinalizeSolutionStep()
            #for node in model_part.Nodes:
            #    node.SetSolutionStepValue(KM.TEMPERATURE, 1.0 )
            

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


    def _GetIOType(self):
        return self.settings["io_settings"]["type"].GetString()

    def __SendControlSignal(self, signal, settings=None):
        data_config = {
            "type"           : "control_signal",
            "control_signal" : signal,
            "settings"       : settings
        }
        self.ExportData(data_config)

    def Finalize(self):
    
        print("control signal : exit")
        self.__SendControlSignal("exit")
        self.vtk_output.ExecuteFinalize()
        super().Finalize()
