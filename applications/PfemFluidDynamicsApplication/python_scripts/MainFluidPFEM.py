from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7

#Activate it to import in the gdb path:
#import sys
#sys.path.append('/home/cpuigbo/kratos')
#x = input("stopped to allow debug: set breakpoints and press enter to continue");
import time as timer
# Import system python
import os

# Import kratos core and applications
import KratosMultiphysics
import KratosMultiphysics.ExternalSolversApplication
import KratosMultiphysics.DelaunayMeshingApplication
import KratosMultiphysics.PfemFluidDynamicsApplication
import KratosMultiphysics.SolidMechanicsApplication

class Solution(object):

    def __init__(self, model, file_parameters = "ProjectParameters.json"):

        self.model=model
        #### TIME MONITORING START ####

        # Time control starts
        print(timer.ctime())
        # Measure process time
        self.t0p = timer.clock()
        # Measure wall time
        self.t0w = timer.time()
        #### TIME MONITORING END ####


        #### PARSING THE PARAMETERS ####

        # Import input
        if file_parameters == "ProjectParameters.json":
            parameter_file = open("ProjectParameters.json",'r')
            self.ProjectParameters = KratosMultiphysics.Parameters(parameter_file.read())
        else:
            self.ProjectParameters = self._import_project_parameters(file_parameters)

        #set echo level
        self.echo_level = self.ProjectParameters["problem_data"]["echo_level"].GetInt()

        # Print solving time
        self.report = False
        if( self.echo_level > 0 ):
            self.report = True

        print(" ")

        # defining the number of threads:
        num_threads = self.ProjectParameters["problem_data"]["threads"].GetInt()
        self.SetParallelSize(num_threads)
        print("::[KPFEM Simulation]:: [OMP USING",num_threads,"THREADS ]")
        #parallel.PrintOMPInfo()


        print(" ")
        print("::[KPFEM Simulation]:: [Time Step:", self.ProjectParameters["problem_data"]["time_step"].GetDouble()," echo:", self.echo_level,"]")

        #### Model_part settings start ####

        # Defining the model_part
        self.main_model_part = self.model.CreateModelPart(self.ProjectParameters["problem_data"]["model_part_name"].GetString())

        self.main_model_part.ProcessInfo.SetValue(KratosMultiphysics.SPACE_DIMENSION, self.ProjectParameters["problem_data"]["dimension"].GetInt())
        self.main_model_part.ProcessInfo.SetValue(KratosMultiphysics.DOMAIN_SIZE, self.ProjectParameters["problem_data"]["dimension"].GetInt())
        self.main_model_part.ProcessInfo.SetValue(KratosMultiphysics.DELTA_TIME, self.ProjectParameters["problem_data"]["time_step"].GetDouble())
        self.main_model_part.ProcessInfo.SetValue(KratosMultiphysics.TIME, self.ProjectParameters["problem_data"]["start_time"].GetDouble())
        if( self.ProjectParameters["problem_data"].Has("gravity_vector") ):
            self.main_model_part.ProcessInfo.SetValue(KratosMultiphysics.GRAVITY_X, self.ProjectParameters["problem_data"]["gravity_vector"][0].GetDouble())
            self.main_model_part.ProcessInfo.SetValue(KratosMultiphysics.GRAVITY_Y, self.ProjectParameters["problem_data"]["gravity_vector"][1].GetDouble())
            self.main_model_part.ProcessInfo.SetValue(KratosMultiphysics.GRAVITY_Z, self.ProjectParameters["problem_data"]["gravity_vector"][2].GetDouble())

        ###TODO replace this "model" for real one once available in kratos core
        self.Model = {self.ProjectParameters["problem_data"]["model_part_name"].GetString() : self.main_model_part}

        #construct the solver (main setting methods are located in the solver_module)
        solver_module = __import__(self.ProjectParameters["solver_settings"]["solver_type"].GetString())
        self.solver = solver_module.CreateSolver(self.main_model_part, self.ProjectParameters["solver_settings"])

        #### Output settings start ####

        self.problem_path = os.getcwd()
        self.problem_name = self.ProjectParameters["problem_data"]["problem_name"].GetString()


    def AddNodalVariablesToModelPart(self):

        # Add variables (always before importing the model part)
        self.solver.AddVariables()

        # Add PfemSolidMechanicsApplication Variables
        import pfem_variables
        pfem_variables.AddVariables(self.main_model_part)


    def Run(self):

        self.Initialize()

        self.RunMainTemporalLoop()

        self.Finalize()


    def Initialize(self):

        # Add variables (always before importing the model part)
        self.AddNodalVariablesToModelPart()

        # Read model_part (note: the buffer_size is set here) (restart is read here)
        self.solver.ImportModelPart()

        # Add dofs (always after importing the model part)
        if((self.main_model_part.ProcessInfo).Has(KratosMultiphysics.IS_RESTARTED)):
            if(self.main_model_part.ProcessInfo[KratosMultiphysics.IS_RESTARTED] == False):
                self.solver.AddDofs()
        else:
            self.solver.AddDofs()

        # Build sub_model_parts or submeshes (rearrange parts for the application of custom processes)
        ## Get the list of the submodel part in the object Model
        for i in range(self.ProjectParameters["solver_settings"]["processes_sub_model_part_list"].size()):
            part_name = self.ProjectParameters["solver_settings"]["processes_sub_model_part_list"][i].GetString()
            if( self.main_model_part.HasSubModelPart(part_name) ):
                self.Model.update({part_name: self.main_model_part.GetSubModelPart(part_name)})

        #### Model_part settings end ####


        #print model_part and properties
        if(self.echo_level>1):
            print("")
            print(self.main_model_part)
            for properties in self.main_model_part.Properties:
                print(properties)

        #### Processes settings start ####

        #obtain the list of the processes to be applied

        import process_handler

        process_parameters = KratosMultiphysics.Parameters("{}")
        process_parameters.AddValue("echo_level", self.ProjectParameters["problem_data"]["echo_level"])
        process_parameters.AddValue("constraints_process_list", self.ProjectParameters["constraints_process_list"])
        process_parameters.AddValue("loads_process_list", self.ProjectParameters["loads_process_list"])
        if( self.ProjectParameters.Has("problem_process_list") ):
            process_parameters.AddValue("problem_process_list", self.ProjectParameters["problem_process_list"])
        if( self.ProjectParameters.Has("output_process_list") ):
            process_parameters.AddValue("output_process_list", self.ProjectParameters["output_process_list"])
        if( self.ProjectParameters.Has("processes_sub_model_part_tree_list") ):
            process_parameters.AddValue("processes_sub_model_part_tree_list",self.ProjectParameters["processes_sub_model_part_tree_list"])
        if( self.ProjectParameters.Has("check_process_list") ):
            process_parameters.AddValue("check_process_list", self.ProjectParameters["check_process_list"])

        self.model_processes = process_handler.ProcessHandler(self.Model, process_parameters)

        self.model_processes.ExecuteInitialize()

        #### processes settings end ####


        # --PLOT GRAPHS OPTIONS START--###############
        #self.problem_path = os.getcwd() #current path
        #plot_active = general_variables.PlotGraphs
        #graph_plot = plot_utils.GraphPlotUtility(model_part, self.problem_path)
        # --PLOT GRAPHS OPTIONS END--#################

        #### START SOLUTION ####

        self.computing_model_part = self.solver.GetComputingModelPart()

        self.graphical_output = self.SetGraphicalOutput()

        ## Sets strategies, builders, linear solvers, schemes and solving info, and fills the buffer
        self.solver.Initialize()
        self.solver.InitializeStrategy()
        self.solver.SetEchoLevel(self.echo_level)



        # Initialize GiD  I/O (gid outputs, file_lists)
        self.GraphicalOutputExecuteInitialize()

        #### Output settings end ####

        # writing a initial state results file
        current_id = 0
        #if(load_restart == False):
        #    if (general_variables.TryToSetTheWeight):
        #        if (general_variables.TryToSetConstantWeight):
        #            conditions.SetConstantWeight( general_variables.TryToSetWeightVertical, general_variables.TryToSetWeightHorizontal);
        #        else:
        #            conditions.SetWeight();

        # set solver info starting parameters
        # solving_info = solving_info_utils.SolvingInfoUtility(model_part, SolverSettings)

        print(" ")
        print("::[KPFEM Simulation]:: Analysis -START- ")

        self.model_processes.ExecuteBeforeSolutionLoop()

        self.GraphicalOutputExecuteBeforeSolutionLoop()

        # Set time settings
        self.step       = self.main_model_part.ProcessInfo[KratosMultiphysics.STEP]
        self.time       = self.main_model_part.ProcessInfo[KratosMultiphysics.TIME]

        self.end_time   = self.ProjectParameters["problem_data"]["end_time"].GetDouble()
        self.delta_time = self.ProjectParameters["problem_data"]["time_step"].GetDouble()


    def RunMainTemporalLoop(self):

        # Solving the problem (time integration)
        while(self.time < self.end_time):

            self.InitializeSolutionStep()
            self.SolveSolutionStep()
            self.FinalizeSolutionStep()


    def InitializeSolutionStep(self):

        self.clock_time = self.StartTimeMeasuring();

        # current time parameters
        # self.main_model_part.ProcessInfo.GetPreviousSolutionStepInfo()[KratosMultiphysics.DELTA_TIME] = self.delta_time
        self.delta_time = self.main_model_part.ProcessInfo[KratosMultiphysics.DELTA_TIME]

        self.time = self.time + self.delta_time
        self.step = self.step + 1

        self.main_model_part.ProcessInfo[KratosMultiphysics.STEP] = self.step
        self.main_model_part.CloneTimeStep(self.time)

        #print(" [STEP:",self.step," TIME:",self.time,"]")

        # processes to be executed at the begining of the solution step
        self.model_processes.ExecuteInitializeSolutionStep()

        self.GraphicalOutputExecuteInitializeSolutionStep()

        # solve time step
        self.solver.InitializeSolutionStep()

        self.StopTimeMeasuring(self.clock_time,"Initialize Step" , self.report);

    def SolveSolutionStep(self):

        self.clock_time = self.StartTimeMeasuring();

        self.solver.Predict()

        self.solver.SolveSolutionStep()

        self.solver.FinalizeSolutionStep()

        self.StopTimeMeasuring(self.clock_time,"Solve Step" , self.report);

    def FinalizeSolutionStep(self):

        self.clock_time = self.StartTimeMeasuring();

        self.GraphicalOutputExecuteFinalizeSolutionStep()

        # processes to be executed at the end of the solution step
        self.model_processes.ExecuteFinalizeSolutionStep()

        #if (self.time>0.0001):
        # processes to be executed before writing the output
        self.model_processes.ExecuteBeforeOutputStep()

        # write output results GiD: (frequency writing is controlled internally)
        self.GraphicalOutputPrintOutput()

        # processes to be executed after witting the output
        self.model_processes.ExecuteAfterOutputStep()

        # Calculate Nodal_Area
        self.CalculateNodalArea()
        
        self.StopTimeMeasuring(self.clock_time,"Finalize Step" , self.report);

    def Finalize(self):

        # Ending the problem (time integration finished)
        self.GraphicalOutputExecuteFinalize()

        self.model_processes.ExecuteFinalize()

        print("::[KPFEM Simulation]:: Analysis -END- ")
        print(" ")

        # Check solving information for any problem
        #~ self.solver.InfoCheck() # InfoCheck not implemented yet.

        #### END SOLUTION ####

        # Measure process time
        tfp = timer.clock()
        # Measure wall time
        tfw = timer.time()

        print("::[KPFEM Simulation]:: [Elapsed Time = %.2f" % (tfw - self.t0w),"seconds] (%.2f" % (tfp - self.t0p),"seconds of cpu/s time)")

        print(timer.ctime())

        # to create a benchmark: add standard benchmark files and decomment next two lines
        # rename the file to: run_test.py
        #from run_test_benchmark_results import *
        #WriteBenchmarkResults(model_part)


    def SetGraphicalOutput(self):
        if( self.ProjectParameters.Has("output_configuration") ):
            from gid_output_process import GiDOutputProcess
            self.output_settings = self.ProjectParameters["output_configuration"]
            return GiDOutputProcess(self.computing_model_part,
                                    self.problem_name,
                                    self.output_settings)
        else:
            return (KratosMultiphysics.Process())

    def GraphicalOutputExecuteInitialize(self):
        self.graphical_output.ExecuteInitialize()

    def GraphicalOutputExecuteBeforeSolutionLoop(self):
        # writing a initial state results file or single file (if no restart)
        if((self.main_model_part.ProcessInfo).Has(KratosMultiphysics.IS_RESTARTED)):
            if(self.main_model_part.ProcessInfo[KratosMultiphysics.IS_RESTARTED] == False):
                self.graphical_output.ExecuteBeforeSolutionLoop()

    def GraphicalOutputExecuteInitializeSolutionStep(self):
        self.graphical_output.ExecuteInitializeSolutionStep()

    def GraphicalOutputExecuteFinalizeSolutionStep(self):
        self.graphical_output.ExecuteFinalizeSolutionStep()

    def GraphicalOutputPrintOutput(self):
        if( self.ProjectParameters.Has("output_configuration") ):
            if(self.graphical_output.IsOutputStep()):
                print("---> Print Output at [STEP:",self.step," TIME:",self.time," DT:",self.delta_time,"]")
                self.graphical_output.PrintOutput()

    def GraphicalOutputExecuteFinalize(self):
        self.graphical_output.ExecuteFinalize()



    def SetParallelSize(self, num_threads):
        parallel = KratosMultiphysics.OpenMPUtils()
        parallel.SetNumThreads(int(num_threads))

    def GetParallelSize(self):
        parallel = KratosMultiphysics.OpenMPUtils()
        return parallel.GetNumThreads()

    def StartTimeMeasuring(self):
        # Measure process time
        time_ip = timer.clock()
        return time_ip

    def StopTimeMeasuring(self, time_ip, process, report):
        # Measure process time
        time_fp = timer.clock()
        if( report ):
            used_time = time_fp - time_ip
            print("::[PFEM Simulation]:: [ %.2f" % round(used_time,2),"s", process," ] ")

    def CalculateNodalArea(self):
        pass

    #### Main internal methods ####

    def _import_project_parameters(self, input_file):
        import input_manager
        self.input_manager = input_manager.InputManager(input_file)
        return self.input_manager.GetProjectParameters()


if __name__ == "__main__":
    model = Model()
    Solution(model).Run()
