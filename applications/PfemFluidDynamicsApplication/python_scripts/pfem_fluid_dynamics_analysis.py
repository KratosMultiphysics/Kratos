from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7

import time as timer
import os
from importlib import import_module

# Import kratos core and applications
import KratosMultiphysics
import KratosMultiphysics.PfemFluidDynamicsApplication as KratosPfemFluid
import KratosMultiphysics.DelaunayMeshingApplication  as KratosDelaunay

from KratosMultiphysics.analysis_stage import AnalysisStage
from KratosMultiphysics.PfemFluidDynamicsApplication import python_solvers_wrapper_pfem_fluid as solver_wrapper

class PfemFluidDynamicsAnalysis(AnalysisStage):
    """The base class for the PfemFluidDynamicsAnalysis
    """
    def __init__(self, model, parameters):
        """The constructor of the AnalysisStage-Object.
        Keyword arguments:
        self -- It signifies an instance of a class.
        model -- The Model to be used
        parameters -- The ProjectParameters used
        """
        self.model = model
        #### TIME MONITORING START ####
        # Time control starts
        self.KratosPrintInfo(timer.ctime())
        # Measure process time
        self.t0p = timer.process_time()
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

        # defining the number of threads:
        num_threads = parameters["problem_data"]["threads"].GetInt()
        self.SetParallelSize(num_threads)
        self.KratosPrintInfo("::[KPFEM Simulation]:: [OMP USING" + str(num_threads) + "THREADS ]")
        #parallel.PrintOMPInfo()

        self.KratosPrintInfo(" ")
        self.KratosPrintInfo("::[KPFEM Simulation]:: [Time Step:" + str(parameters["solver_settings"]["time_stepping"]["time_step"].GetDouble()) + " echo:" +  str(self.echo_level) + "]")

        #### Model_part settings start ####
        super(PfemFluidDynamicsAnalysis,self).__init__(model,parameters)
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

                        #print model_part and properties
        if (self.echo_level>-1):
            for properties in self.main_model_part.Properties:
                self.KratosPrintInfo(properties)

        self.AddPfemVariables()

        if parameters["solver_settings"].Has("constitutive_laws_list"):
            self.constitutive_laws_names     = parameters["solver_settings"]["constitutive_laws_list"]
            self.AddMaterialVariables()
        #else:
        #    self.AddAllMaterialVariables()


    def _CreateSolver(self):
        """Create the solver
        """
        return solver_wrapper.CreateSolverByParameters(self.model, self.project_parameters["solver_settings"],self.project_parameters["problem_data"]["parallel_type"].GetString())

    def Initialize(self):
        """This function initializes the AnalysisStage
        Usage: It is designed to be called ONCE, BEFORE the execution of the solution-loop
        This function has to be implemented in deriving classes!
        """

        # Read model_part from mdpa file
        self._solver.ImportModelPart()

        # Prepare model_part (note: the buffer_size is set here) (restart is read here)
        self._solver.PrepareModelPart()

        # Add dofs (always after importing the model part)
        self._solver.AddDofs()

        #print model_part and properties
        if (self.echo_level>1):
            self.KratosPrintInfo("")
            self.KratosPrintInfo(self.main_model_part)
            for properties in self.main_model_part.Properties:
                self.KratosPrintInfo(properties)

        #### Processes settings start ####

        # obtain the list of the processes to be applied
        from KratosMultiphysics.PfemFluidDynamicsApplication.process_handler import ProcessHandler

        process_parameters = KratosMultiphysics.Parameters("{}")
        process_parameters.AddValue("echo_level", self.project_parameters["problem_data"]["echo_level"])

        if( self.project_parameters.Has("problem_process_list") ):
            process_parameters.AddValue("problem_process_list", self.project_parameters["problem_process_list"])

        self.model_processes = ProcessHandler(self.model, process_parameters)
        self.model_processes.ExecuteInitialize()

        ## here we initialize user-provided processes
        order_processes_initialization = self._GetOrderOfProcessesInitialization()
        self._list_of_processes        = self._CreateProcesses("processes", order_processes_initialization)
        self._list_of_output_processes = self._CreateProcesses("output_processes", order_processes_initialization)
        self._list_of_processes.extend(self._list_of_output_processes)

        for process in self._GetListOfProcesses():
            process.ExecuteInitialize()

        #### processes settings end ####
        #### START SOLUTION ####

        self.computing_model_part = self._solver.GetComputingModelPart()
        self.graphical_output = self.SetGraphicalOutput()
        ## Sets strategies, builders, linear solvers, schemes and solving info, and fills the buffer
        self._solver.Initialize()
        self._solver.InitializeStrategy()
        self._solver.SetEchoLevel(self.echo_level)

        # Initialize GiD  I/O (gid outputs, file_lists)
        self.GraphicalOutputExecuteInitialize()

        self.KratosPrintInfo(" ")
        self.KratosPrintInfo("::[KPFEM Simulation]:: Analysis -START- ")

        self.model_processes.ExecuteBeforeSolutionLoop()

        for process in self._GetListOfProcesses():
            process.ExecuteBeforeSolutionLoop()

        self.GraphicalOutputExecuteBeforeSolutionLoop()

        # write output results GiD: (frequency writing is controlled internally)
        self.GraphicalOutputPrintOutput()

        # Set time settings
        self.step = self.main_model_part.ProcessInfo[KratosMultiphysics.STEP]
        self.time = self.main_model_part.ProcessInfo[KratosMultiphysics.TIME]

        self.end_time   = self.project_parameters["problem_data"]["end_time"].GetDouble()
        self.delta_time = self.project_parameters["solver_settings"]["time_stepping"]["time_step"].GetDouble()


    def InitializeSolutionStep(self):
        """This function performs all the required operations that should be executed
        (for each step) BEFORE solving the solution step.
        """
        self.clock_time = self.StartTimeMeasuring()
        # processes to be executed at the begining of the solution step
        self.model_processes.ExecuteInitializeSolutionStep()

        for process in self._GetListOfProcesses():
            process.ExecuteInitializeSolutionStep()

        self.GraphicalOutputExecuteInitializeSolutionStep()

        # solve time step
        self._solver.InitializeSolutionStep()

        self.StopTimeMeasuring(self.clock_time,"Initialize Step" , self.report)

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
        self.GraphicalOutputPrintOutput()

        # processes to be executed after witting the output
        self.model_processes.ExecuteAfterOutputStep()

        for process in self._GetListOfProcesses():
            process.ExecuteAfterOutputStep()

        self.StopTimeMeasuring(self.clock_time,"Finalize Step" , self.report);

    def Finalize(self):
        """This function finalizes the AnalysisStage
        Usage: It is designed to be called ONCE, AFTER the execution of the solution-loop
        """

        # Ending the problem (time integration finished)
        self.GraphicalOutputExecuteFinalize()
        self.model_processes.ExecuteFinalize()

        for process in self._GetListOfProcesses():
            process.ExecuteFinalize()

        self.KratosPrintInfo("::[KPFEM Simulation]:: Analysis -END- ")
        self.KratosPrintInfo(" ")

        #### END SOLUTION ####

        # Measure process time
        tfp = timer.process_time()
        # Measure wall time
        tfw = timer.time()

        print("::[KPFEM Simulation]:: [Elapsed Time = %.2f" % (tfw - self.t0w),"seconds] (%.2f" % (tfp - self.t0p),"seconds of cpu/s time)")
        self.KratosPrintInfo(timer.ctime())


    def SetGraphicalOutput(self):
        """This function sets the settings for the graphical
        output
        """
        if( self.project_parameters.Has("output_configuration") ):
            from KratosMultiphysics.PfemFluidDynamicsApplication.pfem_fluid_gid_output_process import GiDOutputProcess
            self.output_settings = self.project_parameters["output_configuration"]
            self.post_process_model_part = self.model.CreateModelPart("output_model_part")
            return GiDOutputProcess(self.post_process_model_part,
                                    self.problem_name,
                                    self.output_settings)
        else:
            return (KratosMultiphysics.Process())

    def GraphicalOutputExecuteInitialize(self):
        """This function performs the initialize of the graphical output
        """
        self.graphical_output.ExecuteInitialize()

    def GraphicalOutputExecuteBeforeSolutionLoop(self):
        """This function performs the ExecuteBeforeSolutionLoop
        of the graphical_output
        """
        # writing a initial state results file or single file
        self.graphical_output.ExecuteBeforeSolutionLoop()

    def GraphicalOutputExecuteInitializeSolutionStep(self):
        """This function performs the ExecuteInitializeSolutionStep
        of the graphical_output
        """
        self.graphical_output.ExecuteInitializeSolutionStep()

    def GraphicalOutputExecuteFinalizeSolutionStep(self):
        """This function performs the ExecuteFinalizeSolutionStep
        of the graphical_output
        """
        self.graphical_output.ExecuteFinalizeSolutionStep()

    def GraphicalOutputPrintOutput(self):
        """This function prints the output for this time step
        """
        if( self.project_parameters.Has("output_configuration") ):
            self.post_process_model_part.ProcessInfo[KratosMultiphysics.TIME] = self.main_model_part.ProcessInfo[KratosMultiphysics.TIME]
            if(self.graphical_output.IsOutputStep()):
                time=self.main_model_part.ProcessInfo[KratosMultiphysics.TIME]
                delta_time=self.main_model_part.ProcessInfo[KratosMultiphysics.DELTA_TIME]
                step=self.main_model_part.ProcessInfo[KratosMultiphysics.STEP]
                KratosMultiphysics.PfemFluidDynamicsApplication.PostProcessUtilities().RebuildPostProcessModelPart(self.post_process_model_part, self.main_model_part)
                self.KratosPrintInfo("")
                self.KratosPrintInfo("**********************************************************")
                self.KratosPrintInfo("---> Print Output at [STEP:" + str(step) + " TIME:" + str(time) + " DT:" + str(delta_time) + "]")
                self.KratosPrintInfo("**********************************************************")
                self.KratosPrintInfo("")
                self.graphical_output.PrintOutput()

    def GraphicalOutputExecuteFinalize(self):
        """This function performs the ExecuteFinalize
        of the graphical_output
        """
        self.graphical_output.ExecuteFinalize()

    def SetParallelSize(self, num_threads):
        """This function sets the number of threads
        """
        parallel = KratosMultiphysics.OpenMPUtils()
        parallel.SetNumThreads(int(num_threads))

    def GetParallelSize(self):
        """This function returns the number of threads
        """
        parallel = KratosMultiphysics.OpenMPUtils()
        return parallel.GetNumThreads()

    def StartTimeMeasuring(self):
        """This function starts time calculation
        """
        # Measure process time
        time_ip = timer.process_time()
        return time_ip

    def StopTimeMeasuring(self, time_ip, process, report):
        """This function ends time calculation
        """
        # Measure process time
        time_fp = timer.process_time()
        if report:
            used_time = time_fp - time_ip
            print("::[PFEM Simulation]:: [ %.2f" % round(used_time,2),"s", process," ] ")

    def _GetOrderOfProcessesInitialization(self):
        """This function can be overridden in derived classes if the order of
        initialization for the processes matters
        """
        return ["constraints_process_list",
                "loads_process_list",
                "auxiliar_process_list"]

    def KratosPrintInfo(self, message):
        """This function prints info on screen
        """
        KratosMultiphysics.Logger.Print(message, label="")
        KratosMultiphysics.Logger.Flush()

    def AddMaterialVariables(self):

        if not self.main_model_part.HasNodalSolutionStepVariable(KratosMultiphysics.DENSITY):
            self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.DENSITY)
        if not self.main_model_part.HasNodalSolutionStepVariable(KratosMultiphysics.BULK_MODULUS):
            self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.BULK_MODULUS)
        if not self.main_model_part.HasNodalSolutionStepVariable(KratosMultiphysics.DYNAMIC_VISCOSITY):
            self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.DYNAMIC_VISCOSITY)

        for i in range(self.constitutive_laws_names.size()):
            if (self.constitutive_laws_names[i].GetString()=="FrictionalViscoplastic2DLaw" or self.constitutive_laws_names[i].GetString()=="FrictionalViscoplastic3DLaw"):
                if not self.main_model_part.HasNodalSolutionStepVariable(KratosMultiphysics.INTERNAL_FRICTION_ANGLE):
                    self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.INTERNAL_FRICTION_ANGLE)
                if not self.main_model_part.HasNodalSolutionStepVariable(KratosPfemFluid.COHESION):
                    self.main_model_part.AddNodalSolutionStepVariable(KratosPfemFluid.COHESION)
                if not self.main_model_part.HasNodalSolutionStepVariable(KratosPfemFluid.ADAPTIVE_EXPONENT):
                    self.main_model_part.AddNodalSolutionStepVariable(KratosPfemFluid.ADAPTIVE_EXPONENT)
                if not self.main_model_part.HasNodalSolutionStepVariable(KratosPfemFluid.REGULARIZATION_COEFFICIENT):
                    self.main_model_part.AddNodalSolutionStepVariable(KratosPfemFluid.REGULARIZATION_COEFFICIENT)
            elif (self.constitutive_laws_names[i].GetString()=="Hypoelastic2DLaw" or self.constitutive_laws_names[i].GetString()=="Hypoelastic3DLaw"):
                if not self.main_model_part.HasNodalSolutionStepVariable(KratosMultiphysics.POISSON_RATIO):
                    self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.POISSON_RATIO)
                if not self.main_model_part.HasNodalSolutionStepVariable(KratosMultiphysics.YOUNG_MODULUS):
                    self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.YOUNG_MODULUS)
            elif (self.constitutive_laws_names[i].GetString()=="Bingham2DLaw" or self.constitutive_laws_names[i].GetString()=="Bingham3DLaw"):
                if not self.main_model_part.HasNodalSolutionStepVariable(KratosPfemFluid.FLOW_INDEX):
                    self.main_model_part.AddNodalSolutionStepVariable(KratosPfemFluid.FLOW_INDEX)
                if not self.main_model_part.HasNodalSolutionStepVariable(KratosPfemFluid.YIELD_SHEAR):
                    self.main_model_part.AddNodalSolutionStepVariable(KratosPfemFluid.YIELD_SHEAR)
                if not self.main_model_part.HasNodalSolutionStepVariable(KratosPfemFluid.YIELDED):
                    self.main_model_part.AddNodalSolutionStepVariable(KratosPfemFluid.YIELDED)
                if not self.main_model_part.HasNodalSolutionStepVariable(KratosPfemFluid.ADAPTIVE_EXPONENT):
                    self.main_model_part.AddNodalSolutionStepVariable(KratosPfemFluid.ADAPTIVE_EXPONENT)
            elif (self.constitutive_laws_names[i].GetString()=="PapanastasiouMuIRheology2DLaw" or self.constitutive_laws_names[i].GetString()=="PapanastasiouMuIRheology3DLaw" or
            self.constitutive_laws_names[i].GetString()=="BarkerBercovierMuIRheology2DLaw" or self.constitutive_laws_names[i].GetString()=="BarkerBercovierMuIRheology3DLaw" or
            self.constitutive_laws_names[i].GetString()=="BarkerMuIRheology2DLaw" or self.constitutive_laws_names[i].GetString()=="BarkerMuIRheology3DLaw"):
                if not self.main_model_part.HasNodalSolutionStepVariable(KratosPfemFluid.STATIC_FRICTION):
                    self.main_model_part.AddNodalSolutionStepVariable(KratosPfemFluid.STATIC_FRICTION)
                if not self.main_model_part.HasNodalSolutionStepVariable(KratosPfemFluid.DYNAMIC_FRICTION):
                    self.main_model_part.AddNodalSolutionStepVariable(KratosPfemFluid.DYNAMIC_FRICTION)
                if not self.main_model_part.HasNodalSolutionStepVariable(KratosPfemFluid.INERTIAL_NUMBER_ZERO):
                    self.main_model_part.AddNodalSolutionStepVariable(KratosPfemFluid.INERTIAL_NUMBER_ZERO)
                if not self.main_model_part.HasNodalSolutionStepVariable(KratosPfemFluid.GRAIN_DIAMETER):
                    self.main_model_part.AddNodalSolutionStepVariable(KratosPfemFluid.GRAIN_DIAMETER)
                if not self.main_model_part.HasNodalSolutionStepVariable(KratosPfemFluid.GRAIN_DENSITY):
                    self.main_model_part.AddNodalSolutionStepVariable(KratosPfemFluid.GRAIN_DENSITY)
                if not self.main_model_part.HasNodalSolutionStepVariable(KratosPfemFluid.REGULARIZATION_COEFFICIENT):
                    self.main_model_part.AddNodalSolutionStepVariable(KratosPfemFluid.REGULARIZATION_COEFFICIENT)
                if (self.constitutive_laws_names[i].GetString()=="BarkerBercovierMuIRheology2DLaw" or self.constitutive_laws_names[i].GetString()=="BarkerBercovierMuIRheology3DLaw" or
                self.constitutive_laws_names[i].GetString()=="BarkerMuIRheology2DLaw" or self.constitutive_laws_names[i].GetString()=="BarkerMuIRheology3DLaw"):
                    if not self.main_model_part.HasNodalSolutionStepVariable(KratosPfemFluid.INFINITE_FRICTION):
                        self.main_model_part.AddNodalSolutionStepVariable(KratosPfemFluid.INFINITE_FRICTION)
                    if not self.main_model_part.HasNodalSolutionStepVariable(KratosPfemFluid.INERTIAL_NUMBER_ONE):
                        self.main_model_part.AddNodalSolutionStepVariable(KratosPfemFluid.INERTIAL_NUMBER_ONE)
                    if not self.main_model_part.HasNodalSolutionStepVariable(KratosPfemFluid.ALPHA_PARAMETER):
                        self.main_model_part.AddNodalSolutionStepVariable(KratosPfemFluid.ALPHA_PARAMETER)
            elif (self.constitutive_laws_names[i].GetString()!="None" and self.constitutive_laws_names[i].GetString()!="Newtonian2DLaw" and self.constitutive_laws_names[i].GetString()!="Newtonian3DLaw"):
                print("ERROR: THE CONSTITUTIVE LAW PROVIDED FOR THIS SUBMODEL PART IS NOT IN THE PFEM FLUID DATABASE")


    def AddAllMaterialVariables(self):
        print("ATTENTION! YOU ARE ADDING ALL MATERIAL VARIABLES, PLEASE UPDATE YOUR PROJECTPARAMETERS.JSON")

        if not self.main_model_part.HasNodalSolutionStepVariable(KratosMultiphysics.DENSITY):
            self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.DENSITY)
        if not self.main_model_part.HasNodalSolutionStepVariable(KratosMultiphysics.BULK_MODULUS):
            self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.BULK_MODULUS)
        if not self.main_model_part.HasNodalSolutionStepVariable(KratosMultiphysics.DYNAMIC_VISCOSITY):
            self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.DYNAMIC_VISCOSITY)
        if not self.main_model_part.HasNodalSolutionStepVariable(KratosMultiphysics.INTERNAL_FRICTION_ANGLE):
            self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.INTERNAL_FRICTION_ANGLE)
        if not self.main_model_part.HasNodalSolutionStepVariable(KratosPfemFluid.COHESION):
            self.main_model_part.AddNodalSolutionStepVariable(KratosPfemFluid.COHESION)
        if not self.main_model_part.HasNodalSolutionStepVariable(KratosPfemFluid.ADAPTIVE_EXPONENT):
            self.main_model_part.AddNodalSolutionStepVariable(KratosPfemFluid.ADAPTIVE_EXPONENT)
        if not self.main_model_part.HasNodalSolutionStepVariable(KratosPfemFluid.REGULARIZATION_COEFFICIENT):
            self.main_model_part.AddNodalSolutionStepVariable(KratosPfemFluid.REGULARIZATION_COEFFICIENT)
        if not self.main_model_part.HasNodalSolutionStepVariable(KratosMultiphysics.POISSON_RATIO):
            self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.POISSON_RATIO)
        if not self.main_model_part.HasNodalSolutionStepVariable(KratosMultiphysics.YOUNG_MODULUS):
            self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.YOUNG_MODULUS)
        if not self.main_model_part.HasNodalSolutionStepVariable(KratosPfemFluid.FLOW_INDEX):
            self.main_model_part.AddNodalSolutionStepVariable(KratosPfemFluid.FLOW_INDEX)
        if not self.main_model_part.HasNodalSolutionStepVariable(KratosPfemFluid.YIELD_SHEAR):
            self.main_model_part.AddNodalSolutionStepVariable(KratosPfemFluid.YIELD_SHEAR)
        if not self.main_model_part.HasNodalSolutionStepVariable(KratosPfemFluid.YIELDED):
            self.main_model_part.AddNodalSolutionStepVariable(KratosPfemFluid.YIELDED)
        if not self.main_model_part.HasNodalSolutionStepVariable(KratosPfemFluid.STATIC_FRICTION):
            self.main_model_part.AddNodalSolutionStepVariable(KratosPfemFluid.STATIC_FRICTION)
        if not self.main_model_part.HasNodalSolutionStepVariable(KratosPfemFluid.DYNAMIC_FRICTION):
            self.main_model_part.AddNodalSolutionStepVariable(KratosPfemFluid.DYNAMIC_FRICTION)
        if not self.main_model_part.HasNodalSolutionStepVariable(KratosPfemFluid.INERTIAL_NUMBER_ZERO):
            self.main_model_part.AddNodalSolutionStepVariable(KratosPfemFluid.INERTIAL_NUMBER_ZERO)
        if not self.main_model_part.HasNodalSolutionStepVariable(KratosPfemFluid.GRAIN_DIAMETER):
            self.main_model_part.AddNodalSolutionStepVariable(KratosPfemFluid.GRAIN_DIAMETER)
        if not self.main_model_part.HasNodalSolutionStepVariable(KratosPfemFluid.GRAIN_DENSITY):
            self.main_model_part.AddNodalSolutionStepVariable(KratosPfemFluid.GRAIN_DENSITY)
        if not self.main_model_part.HasNodalSolutionStepVariable(KratosPfemFluid.INFINITE_FRICTION):
            self.main_model_part.AddNodalSolutionStepVariable(KratosPfemFluid.INFINITE_FRICTION)
        if not self.main_model_part.HasNodalSolutionStepVariable(KratosPfemFluid.INERTIAL_NUMBER_ONE):
            self.main_model_part.AddNodalSolutionStepVariable(KratosPfemFluid.INERTIAL_NUMBER_ONE)
        if not self.main_model_part.HasNodalSolutionStepVariable(KratosPfemFluid.ALPHA_PARAMETER):
            self.main_model_part.AddNodalSolutionStepVariable(KratosPfemFluid.ALPHA_PARAMETER)


    def AddPfemVariables(self):

        if not self.main_model_part.HasNodalSolutionStepVariable(KratosMultiphysics.MESH_VELOCITY):
            self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.MESH_VELOCITY)

        if not self.main_model_part.HasNodalSolutionStepVariable(KratosMultiphysics.BODY_FORCE):
            self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.BODY_FORCE)

        if not self.main_model_part.HasNodalSolutionStepVariable(KratosMultiphysics.NODAL_MASS):
            self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.NODAL_MASS)

        if not self.main_model_part.HasNodalSolutionStepVariable(KratosMultiphysics.REACTION):
            self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.REACTION)

        if not self.main_model_part.HasNodalSolutionStepVariable(KratosMultiphysics.NORMAL):
            self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.NORMAL)

        if not self.main_model_part.HasNodalSolutionStepVariable(KratosMultiphysics.VOLUME_ACCELERATION):
            self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.VOLUME_ACCELERATION)

        if not self.main_model_part.HasNodalSolutionStepVariable(KratosPfemFluid.FREESURFACE):
            self.main_model_part.AddNodalSolutionStepVariable(KratosPfemFluid.FREESURFACE)

        if not self.main_model_part.HasNodalSolutionStepVariable(KratosPfemFluid.PREVIOUS_FREESURFACE):
            self.main_model_part.AddNodalSolutionStepVariable(KratosPfemFluid.PREVIOUS_FREESURFACE)

        if not self.main_model_part.HasNodalSolutionStepVariable(KratosPfemFluid.PRESSURE_VELOCITY):
            self.main_model_part.AddNodalSolutionStepVariable(KratosPfemFluid.PRESSURE_VELOCITY)

        if not self.main_model_part.HasNodalSolutionStepVariable(KratosPfemFluid.PRESSURE_ACCELERATION):
             self.main_model_part.AddNodalSolutionStepVariable(KratosPfemFluid.PRESSURE_ACCELERATION)

        if not self.main_model_part.HasNodalSolutionStepVariable(KratosPfemFluid.ISOLATED_NODE):
            self.main_model_part.AddNodalSolutionStepVariable(KratosPfemFluid.ISOLATED_NODE)

        if not self.main_model_part.HasNodalSolutionStepVariable(KratosMultiphysics.NODAL_H):
            self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.NODAL_H)

        if not self.main_model_part.HasNodalSolutionStepVariable(KratosDelaunay.SHRINK_FACTOR):
            self.main_model_part.AddNodalSolutionStepVariable(KratosDelaunay.SHRINK_FACTOR)

        if not self.main_model_part.HasNodalSolutionStepVariable(KratosDelaunay.PROPERTY_ID):
            self.main_model_part.AddNodalSolutionStepVariable(KratosDelaunay.PROPERTY_ID)



if __name__ == "__main__":
    parameter_file_name = "ProjectParameters.json"
    with open(parameter_file_name,'r') as parameter_file:
        parameters = Kratos.Parameters(parameter_file.read())
    model = KratosMultiphysics.Model()
    simulation = PfemFluidDynamicsAnalysis(model,parameters)
    simulation.Run()
