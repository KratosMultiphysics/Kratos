
import KratosMultiphysics
import time as timer
import KratosMultiphysics.PfemFluidDynamicsApplication.pfem_fluid_dynamics_analysis as PfemFluidDynamicsAnalysis
import os
from importlib import import_module

def Wait():
    input("Alejandro -> Press Something")

def KratosPrintInfo(message):
    KratosMultiphysics.Logger.Print(message, label="")
    KratosMultiphysics.Logger.Flush()

#============================================================================================================================
class MainPFEM_for_coupling_solution(PfemFluidDynamicsAnalysis.PfemFluidDynamicsAnalysis):
    """
    The derived class for the PfemFluidDynamicsAnalysis
    """
#============================================================================================================================
#============================================================================================================================
    def __init__(self, model, FEM_model_part, parameters):
        """
        The constructor of the MainPFEM_for_coupling_solution-Object.
        Keyword arguments:
        self -- It signifies an instance of a class.
        model -- The Model to be used
        parameters -- The ProjectParameters used
        """

        # We change the name of the inputs in order to be separated from FEMDEM
        problem_name = parameters["problem_data"]["problem_name"].GetString()
        parameters["problem_data"]["problem_name"].SetString("PFEM" + problem_name)
        parameters["solver_settings"]["model_import_settings"]["input_filename"].SetString("PFEM" + problem_name)

        folders = problem_name.split("/")
        if len(folders) > 1:
            parameters["solver_settings"]["model_import_settings"]["input_filename"].SetString(problem_name)

        self.FEM_model_part = FEM_model_part

        self.model = model
        #### TIME MONITORING START ####
        # Time control starts
        self.KratosPrintInfo(timer.ctime())
        # Measure process time
        self.t0p = timer.clock()
        # Measure wall time
        self.t0w = timer.time()
        #### TIME MONITORING END ####

        #### PARSING THE PARAMETERS ####

        #set echo level
        self.echo_level = parameters["problem_data"]["echo_level"].GetInt()

        # Print solving time
        self.report = False
        if( self.echo_level > 0 ):
            self.report = True

        self.KratosPrintInfo(" ")
        self.KratosPrintInfo("::[KPFEM Simulation]:: [Time Step:" + str(parameters["solver_settings"]["time_stepping"]["time_step"].GetDouble()) + " echo:" +  str(self.echo_level) + "]")

        #### Model_part settings start ####
        self.model = model
        self.project_parameters = parameters

        ## Get echo level and parallel type
        self.echo_level = self.project_parameters["problem_data"]["echo_level"].GetInt()
        self.parallel_type = self.project_parameters["problem_data"]["parallel_type"].GetString()
        is_distributed_run = KratosMultiphysics.IsDistributedRun()

        if self.parallel_type == "OpenMP" and is_distributed_run:
            KratosMultiphysics.Logger.PrintWarning("Parallel Type", '"OpenMP" is specified as "parallel_type", but Kratos is running distributed!')
        if self.parallel_type == "MPI" and not is_distributed_run:
            KratosMultiphysics.Logger.PrintWarning("Parallel Type", '"MPI" is specified as "parallel_type", but Kratos is not running distributed!')

        self._GetSolver().AddVariables() # this creates the solver and adds the variables
        # Defining the model_part
        self.main_model_part = self.model.GetModelPart(parameters["solver_settings"]["model_part_name"].GetString())

        self.main_model_part.ProcessInfo.SetValue(KratosMultiphysics.SPACE_DIMENSION, parameters["solver_settings"]["domain_size"].GetInt())
        self.main_model_part.ProcessInfo.SetValue(KratosMultiphysics.DOMAIN_SIZE, parameters["solver_settings"]["domain_size"].GetInt())
        self.main_model_part.ProcessInfo.SetValue(KratosMultiphysics.DELTA_TIME, parameters["solver_settings"]["time_stepping"]["time_step"].GetDouble())
        self.main_model_part.ProcessInfo.SetValue(KratosMultiphysics.TIME, parameters["problem_data"]["start_time"].GetDouble())
        if parameters["problem_data"].Has("gravity_vector"):
             self.main_model_part.ProcessInfo.SetValue(KratosMultiphysics.GRAVITY_X, parameters["problem_data"]["gravity_vector"][0].GetDouble())
             self.main_model_part.ProcessInfo.SetValue(KratosMultiphysics.GRAVITY_Y, parameters["problem_data"]["gravity_vector"][1].GetDouble())
             self.main_model_part.ProcessInfo.SetValue(KratosMultiphysics.GRAVITY_Z, parameters["problem_data"]["gravity_vector"][2].GetDouble())

        self.problem_path = os.getcwd()
        self.problem_name = parameters["problem_data"]["problem_name"].GetString()

#============================================================================================================================
    def SetCustomGraphicalOutput(self, custom_parameters):
        """This function sets the settings for the graphical
        output
        """
        if custom_parameters.Has("output_configuration"):
            from KratosMultiphysics.PfemFluidDynamicsApplication.pfem_fluid_gid_output_process import GiDOutputProcess
            self.output_settings = custom_parameters["output_configuration"]
            self.output_settings["result_file_configuration"].RemoveValue("nodal_results")
            self.output_settings["result_file_configuration"].AddValue("nodal_results", self.project_parameters["output_configuration"]["result_file_configuration"]["nodal_results"])
            self.output_settings["result_file_configuration"].RemoveValue("gauss_point_results")
            self.output_settings["result_file_configuration"].AddValue("gauss_point_results", self.project_parameters["output_configuration"]["result_file_configuration"]["gauss_point_results"])
            return GiDOutputProcess(self.post_process_model_part,
                                    self.problem_name,
                                    self.output_settings)
        else:
            return (KratosMultiphysics.Process())

#============================================================================================================================
    def _CreateSolver(self):
        """Create the solver
        """
        python_module_name = "KratosMultiphysics.FemToDemApplication"
        full_module_name = python_module_name + "." + "pfem_fluid_solver_for_coupling"
        solver_module = import_module(full_module_name)
        solver = solver_module.CreateSolver(self.model, self.FEM_model_part, self.project_parameters["solver_settings"])
        return solver
        
#============================================================================================================================

    def FinalizeSolutionStep(self):
        """This function performs all the required operations that should be executed
        (for each step) AFTER solving the solution step.
        """
        self.clock_time = self.StartTimeMeasuring();
        self._GetSolver().FinalizeSolutionStep()
        self.GraphicalOutputExecuteFinalizeSolutionStep()

        # processes to be executed at the end of the solution step
        self.model_processes.ExecuteFinalizeSolutionStep()

        for process in self._GetListOfProcesses():
            process.ExecuteFinalizeSolutionStep()
        self.model_processes.ExecuteBeforeOutputStep()

        for process in self._GetListOfProcesses():
            process.ExecuteBeforeOutputStep()

        # write output results GiD: (frequency writing is controlled internally)
        # self.GraphicalOutputPrintOutput()

        # processes to be executed after witting the output
        self.model_processes.ExecuteAfterOutputStep()

        for process in self._GetListOfProcesses():
            process.ExecuteAfterOutputStep()

        self.StopTimeMeasuring(self.clock_time,"Finalize Step" , self.report);