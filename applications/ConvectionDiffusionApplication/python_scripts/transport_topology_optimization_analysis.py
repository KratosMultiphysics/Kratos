# Importing the Kratos Library
import KratosMultiphysics


# Import applications
import KratosMultiphysics.ConvectionDiffusionApplication as KratosCD

from KratosMultiphysics.ConvectionDiffusionApplication.convection_diffusion_analysis import ConvectionDiffusionAnalysis
from KratosMultiphysics.analysis_stage import AnalysisStage
from KratosMultiphysics.ConvectionDiffusionApplication import transport_topology_optimization_solver

class TransportTopologyOptimizationAnalysis(ConvectionDiffusionAnalysis):
    def __init__(self,model,parameters, constant_conv_vel=True):
        self.topology_optimization_stage = 0
        self.topology_optimization_stage_str = "INIT"
        super().__init__(model,parameters)
        self._DefineTransportProperties() 
        self._DefineConvectiveVelocity(constant_conv_vel)
        
    def _DefineTransportProperties(self, decay = 0.0, convection_coefficient = 1.0):
        self._GetTransportSolver()._DefineProperties(decay, convection_coefficient)
        self._GetAdjointTransportSolver()._DefineProperties(decay, convection_coefficient)

    def _DefineConvectiveVelocity(self, const_convective_vel, velocity=[0,0,0]):
        self._GetTransportSolver()._DefineConvectionVelocity(const_convective_vel, velocity)
        self._GetAdjointTransportSolver()._DefineConvectionVelocity(const_convective_vel, velocity)

    def _CheckMaterialProperties(self):
        self._GetSolver()._CheckMaterialProperties()

    def _CreateSolver(self, isAdjointSolver = False):
        """
        This method creates a solver
        isAdjointSolver == False --> T_solver
        isAdjointSolver == True  --> ADJ_T_solver
        """
        return transport_topology_optimization_solver.CreateSolver(self.model, self.project_parameters, isAdjointSolver)

    def RunSolutionLoop(self):
        # self.CheckProcesses()
        self._SolveTransportProblem()
        self._SolveAdjointTransportProblem()

    def _RunStageSolutionLoop(self, problem_stage):
        self._SetTopologyOptimizationStage(problem_stage)
        # self._CheckMaterialProperties()
        self._PrepareSettings()
        # input("material prop")
        super().RunSolutionLoop()

    def _SolveTransportProblem(self):
        self._RunStageSolutionLoop(problem_stage=1)

    def _SolveAdjointTransportProblem(self):
        self._RunStageSolutionLoop(problem_stage=2)

    def CheckProcesses(self):
        """here the boundary conditions is applied, by calling the InitializeSolutionStep function of the processes"""

        for process in self._GetListOfProcesses():
            name = type(process).__name__
            if hasattr(process, "variable"):
                variable = process.variable
            else:
                variable = "Not defined"
            if hasattr(process, "value"):
                value = process.value
            elif hasattr(process, "function_string"):
                value = process.function_string
            else:
                value = "Not defined"
            print("PROCESS:", name, ", VAR:", variable, ", VALUE:", value)
            input()
    
    def _GetSolver(self, force_adjoint = False):
        """
        This method returns solver in use for the current topology optimization phase
        """
        if not hasattr(self, 'T_solver'):
                self.T_solver = self._CreateSolver()
        if not hasattr(self, 'ADJ_T_solver'):
                self.ADJ_T_solver = self._CreateSolver(True)
        self._solver = self._GetTopologyOptimizationStageSolver(force_adjoint)
        return self._solver
    
    def _GetTransportSolver(self):
        """
        This method returns the Transport Solver
        """
        if not hasattr(self, 'T_solver'):
                self.T_solver = self._CreateSolver(False)
        return self.T_solver
    
    def _GetAdjointTransportSolver(self):
        """
        This method returns the Adjoint Transport Solver
        """
        if not hasattr(self, 'ADJ_T_solver'):
                self.ADJ_T_solver = self._CreateSolver(True)
        return self.ADJ_T_solver
    
    def _GetTopologyOptimizationStageSolver(self, force_adjont = False):
        """
        This methods returns the current topology optimization stage solver
        iF force_adjoint --> return ADJ_T_solver
        If topology_optimization_stage != 2 (<=> EVERYTHING BUT NOT ADJ STAGE) --> return T_solver
        If topology_optimization_stage == 2 (<=>  ADJ STAGE) --> return ADJ_T_solver
        """
        if (self.IsAdjointTransportStage() or (force_adjont)): # ADJ
            return self.ADJ_T_solver
        else: # NS
            return self.T_solver
        
    def IsInitializeStage(self):
        """
        This method returns TRUE the topology optimization stage is 0 <=> INIT
        """ 
        return self.CheckOptimizationStage(0)
    
    def IsTransportStage(self):
        """
        This method returns TRUE the topology optimization stage is 1 <=> NS
        """
        return self.CheckOptimizationStage(1)

    def IsAdjointTransportStage(self):
        """
        This method returns TRUE the topology optimization stage is 2 <=> ADJ
        """
        return self.CheckOptimizationStage(2)
    
    def GetOptimizationStage(self):
        """
        This method returns the topology optimization stage
        """
        return self.topology_optimization_stage

    def CheckOptimizationStage(self, check):
        """
        This method returns TRUE the topology optimization stage is 'check'
        """ 
        return (self.GetOptimizationStage() == check)
    
    def Initialize(self):
        """This method initializes the FluidTopologyOptimizationAnalysis
        Usage: It is designed to be called ONCE, BEFORE the execution of the solution-loop
        """
        # Modelers:
        self._CreateModelers()
        self._ModelersSetupGeometryModel()
        self._ModelersPrepareGeometryModel()
        self._ModelersSetupModelPart()
        # Prepare Solvers : ImportModelPart -> PrepareModelPart -> AddDofs)
        self.PrepareSolvers()
        # Set the Functionals Weigths for the Optimization
        self._SetFunctionalWeights()
        # Modify Initial
        self.ModifyInitialProperties()
        self.ModifyInitialGeometry()
        # Initialize user-provided processes
        self.__CreateListOfProcesses() # has to be done after importing and preparing the ModelPart
        for process in self._GetListOfProcesses():
            process.ExecuteInitialize()
        # Initialize Solvers
        self.InitializeSolvers()
        # Checks
        self.Check()
        self.ModifyAfterSolverInitialize()
        # Set Topology Optimization Stage: Initialize
        self._SetTopologyOptimizationStage(0)
        for process in self._GetListOfProcesses():
            process.ExecuteBeforeSolutionLoop()
        ## Stepping and time settings
        if self._GetComputingModelPart().ProcessInfo[KratosMultiphysics.IS_RESTARTED]:
            self.time = self._GetComputingModelPart().ProcessInfo[KratosMultiphysics.TIME]
        else:
            self.time = self.project_parameters["problem_data"]["start_time"].GetDouble()
            self._GetComputingModelPart().ProcessInfo[KratosMultiphysics.TIME] = self.time
        self.start_time = self.time
        self.end_time = self.project_parameters["problem_data"]["end_time"].GetDouble()
        ## If the echo level is high enough, print the complete list of settings used to run the simulation
        if self.echo_level > 1:
            with open("ProjectParametersOutput.json", 'w') as parameter_output_file:
                parameter_output_file.write(self.project_parameters.PrettyPrintJsonString())
        # Print Start Analysis
        KratosMultiphysics.Logger.PrintInfo("\n" + self._GetSimulationName(), "Analysis -START- ")

    def PrepareSolvers(self):
        """This method prepares the NS and ADJ_NS Solvers in the AnalysisStage 
        Usage: It is designed to be called ONCE, BEFORE the execution of the solution-loop
        Prepare Solvers : ImportModelPart -> PrepareModelPart -> AddDofs
        """
        self.PrepareTransportSolver()
        self.PrepareAdjointTransportSolver()
    
    def PrepareTransportSolver(self):
        """This method prepares the Navier-Stokes primal problem Solver in the AnalysisStage 
        Usage: It is designed to be called ONCE, BEFORE the execution of the solution-loop
        Prepare Solver : ImportModelPart -> PrepareModelPart -> AddDofs
        """
        self._GetTransportSolver().ImportModelPart()
        self._GetTransportSolver().PrepareModelPart()
        self._GetTransportSolver().AddDofs()
    
    def PrepareAdjointTransportSolver(self):
        """This method prepares the Adjoint Navier-Stokes Solver in the AnalysisStage 
        Usage: It is designed to be called ONCE, BEFORE the execution of the solution-loop
        Prepare Solver : ImportModelPart -> PrepareModelPart -> AddDofs
        """
        # Modelers:
        self._GetAdjointTransportSolver().ImportModelPart()
        self._GetAdjointTransportSolver().PrepareModelPart()
        self._GetAdjointTransportSolver().AddDofs()

    def InitializeSolvers(self):
        """This method initializes the NS and ADJ_NS Solvers in the AnalysisStage 
        Usage: It is designed to be called ONCE, BEFORE the execution of the solution-loop
        """
        self.InitializeTransportSolver()
        self.InitializeAdjointTransportSolver()

    def InitializeTransportSolver(self):
        """This method initializes the NS Solver in the AnalysisStage 
        Usage: It is designed to be called ONCE, BEFORE the execution of the solution-loop
        """
        self._GetTransportSolver().Initialize()

    def InitializeAdjointTransportSolver(self):
        """This method initializes the ADJ_NS Solver in the AnalysisStage 
        Usage: It is designed to be called ONCE, BEFORE the execution of the solution-loop
        """
        self._GetAdjointTransportSolver().Initialize()

    def _SetFunctionalWeights(self):
        dumb = 0
        return dumb
    
    def __CreateListOfProcesses(self):
        """This method creates the processes and the output-processes
        """
        order_processes_initialization = self._GetOrderOfProcessesInitialization()
        self._list_of_processes        = self._CreateProcesses("processes", order_processes_initialization)
        deprecated_output_processes    = self._CheckDeprecatedOutputProcesses(self._list_of_processes)
        order_processes_initialization = self._GetOrderOfOutputProcessesInitialization()
        self._list_of_output_processes = self._CreateProcesses("output_processes", order_processes_initialization)
        self._list_of_processes.extend(self._list_of_output_processes) # Adding the output processes to the regular processes
        self._list_of_output_processes.extend(deprecated_output_processes)

    def _SetTopologyOptimizationStage(self, problem_stage):
        """
        This method set the topology optimization problem phase:
        0: initialization, 1: NS solution, 2: ADJ NS solution, 3: Optimization Update
        """
        self.topology_optimization_stage = problem_stage 
        if self.CheckOptimizationStage(0):
            self.topology_optimization_stage_str = "INI"
        elif self.CheckOptimizationStage(1):
            self.topology_optimization_stage_str = " T "
        elif self.CheckOptimizationStage(2):
            self.topology_optimization_stage_str = "A_T"
        elif self.CheckOptimizationStage(3):
            self.topology_optimization_stage_str = "OPT"
        else:
            self.topology_optimization_stage_str = "ERROR"
        self._GetSolver()._SetTopologyOptimizationStage(self.topology_optimization_stage)

    def _GetComputingModelPart(self):
        """
        This method returns the computing model part for the current solver
        """    
        return self._GetSolver().GetComputingModelPart()
    
    def KeepAdvancingSolutionLoop(self):
        """This method specifies the stopping criteria for breaking the solution loop
        It can be overridden by derived classes
        """
        if self.IsTransportStage(): # T
            return self.time < self.end_time
        elif self.IsAdjointTransportStage(): #ADJ
            return self.time > self.start_time
        else:
            print("--|" + self.topology_optimization_stage_str + "| Time check outside T or ADJ_T solution")
            return False
        
    def _AdvanceTime(self):
        """
        This methods hadls the time for the Topology Optimization problems, currently its purpos it's just to make thinks work
        """
        time = super()._AdvanceTime()
        if self.IsAdjointTransportStage():
            time = 0.0
        return time
    
    def _PrepareSettings(self):
        if (self.IsTransportStage()):
            # self._GetTransportSolver().settings["convection_diffusion_variables"]["unknown_variable"].SetString("TEMPERATURE")
            # self._GetTransportSolver().settings["convection_diffusion_variables"]["volume_source_variable"].SetString("HEAT_FLUX")
            # self._GetTransportSolver().settings["convection_diffusion_variables"]["surface_source_variable"].SetString("FACE_HEAT_FLUX")
            convention_diffusion_settings = self._GetComputingModelPart().ProcessInfo.GetValue(KratosMultiphysics.CONVECTION_DIFFUSION_SETTINGS)
            convention_diffusion_settings.SetUnknownVariable(KratosMultiphysics.KratosGlobals.GetVariable("TEMPERATURE"))
            convention_diffusion_settings.SetVolumeSourceVariable(KratosMultiphysics.KratosGlobals.GetVariable("HEAT_FLUX"))
            convention_diffusion_settings.SetSurfaceSourceVariable(KratosMultiphysics.KratosGlobals.GetVariable("FACE_HEAT_FLUX"))
        elif (self.IsAdjointTransportStage()):
            # self._GetAdjointTransportSolver().settings["convection_diffusion_variables"]["unknown_variable"].SetString("TEMPERATURE_ADJ")
            # self._GetAdjointTransportSolver().settings["convection_diffusion_variables"]["volume_source_variable"].SetString("HEAT_FLUX_ADJ")
            # self._GetAdjointTransportSolver().settings["convection_diffusion_variables"]["surface_source_variable"].SetString("FACE_HEAT_FLUX_ADJ")
            convention_diffusion_settings = self._GetComputingModelPart().ProcessInfo.GetValue(KratosMultiphysics.CONVECTION_DIFFUSION_SETTINGS)
            convention_diffusion_settings.SetUnknownVariable(KratosMultiphysics.KratosGlobals.GetVariable("TEMPERATURE_ADJ"))
            convention_diffusion_settings.SetVolumeSourceVariable(KratosMultiphysics.KratosGlobals.GetVariable("HEAT_FLUX_ADJ"))
            convention_diffusion_settings.SetSurfaceSourceVariable(KratosMultiphysics.KratosGlobals.GetVariable("FACE_HEAT_FLUX_ADJ"))
        else:
            pass