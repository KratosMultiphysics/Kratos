# Importing the Kratos Library
import KratosMultiphysics

import numpy as np #import the numpy library
import scipy as sp #import the scipy library


# Import applications
import KratosMultiphysics.ConvectionDiffusionApplication as KratosCD

from KratosMultiphysics.ConvectionDiffusionApplication.convection_diffusion_analysis import ConvectionDiffusionAnalysis
from KratosMultiphysics.analysis_stage import AnalysisStage
from KratosMultiphysics.ConvectionDiffusionApplication import transport_topology_optimization_solver

class TransportTopologyOptimizationAnalysisTest(ConvectionDiffusionAnalysis):
    def __init__(self,model,parameters):
        self.topology_optimization_stage = 0
        self.topology_optimization_stage_str = "INIT"
        super().__init__(model,parameters)
        # self._DefineTransportProperties() 
        
    # def _DefineTransportProperties(self, decay = 0.0, convection_coefficient = 1.0, const_velocity=False, vel=[0,0,0]):
    #     self._GetPhysicsSolver()._DefineProperties(decay, convection_coefficient, const_velocity, vel)
    #     self._GetAdjointSolver()._DefineProperties(decay, convection_coefficient, const_velocity, vel)

    def _CheckMaterialProperties(self):
        print("--|CHECK| Check Transport Properties")
        self._GetSolver()._CheckMaterialProperties()

    def _CreateSolver(self, isAdjointSolver = False):
        """
        This method creates a solver
        isAdjointSolver == False --> physics_solver
        isAdjointSolver == True  --> adjoint_solver
        """
        return transport_topology_optimization_solver.CreateSolver(self.model, self.project_parameters, isAdjointSolver)
    
    def _InitializePhysicsParameters(self, conductivity_parameters=[1e-4, 1e-2, 1.0], decay_parameters=[0.0, 0.0, 1.0], convection_coefficient_parameters=[0.0, 0.0, 1.0]):
        print("--|" + self.topology_optimization_stage_str + "| INITIALIZE PHYSICS PARAMETERS")
        self._InitializeConductivity(conductivity_parameters)
        self._InitializeDecay(decay_parameters)
        self._InitializeConvectionCoefficient(convection_coefficient_parameters)

    def _InitializeConductivity(self, conductivity_parameters_values):
        print("--|" + self.topology_optimization_stage_str + "| INITIALIZE CONDUCTIVITY")
        self.conductivity_parameters = conductivity_parameters_values
        self._ResetConductivity()
        self._UpdateConductivityVariable()

    def _InitializeDecay(self, decay_parameters_values):
        print("--|" + self.topology_optimization_stage_str + "| INITIALIZE DECAY")
        self.decay_parameters = decay_parameters_values
        self._ResetDecay()
        self._UpdateDecayVariable()

    def _InitializeConvectionCoefficient(self, convection_coefficient_parameters_values):
        print("--|" + self.topology_optimization_stage_str + "| INITIALIZE CONVECTION COEFFICIENT")
        self.convection_coefficient_parameters = convection_coefficient_parameters_values
        self._ResetConvectionCoefficient()
        self._SetConvectiveVelocity()
        self._UpdateConvectionCoefficientVariable()

    def ResetPhysicsParameters(self):
        self._ResetConductivity()
        self._ResetDecay()
        self._ResetConvectionCoefficient()

    def _ResetConductivity(self):
        self.conductivity = 0.0

    def _ResetDecay(self):
        self.decay = 0.0

    def _ResetConvectionCoefficient(self):
        self.convection_coefficient = 0.0

    def _SetConvectiveVelocity(self, constant_velocity=False, vel=[0,0,0]):
        self.is_constant_velocity = constant_velocity
        self.constant_velocity = vel
        self._GetSolver()._SetConvectiveVelocity(self.is_constant_velocity, self.constant_velocity)

    def RunSolutionLoop(self):
        # self.CheckProcesses()
        self.first_iteration = True
        self._InitializePhysicsParameters()
        self._InitializeTopologyOptimizationStepPhysicsSolution()
        self._CheckMaterialProperties()
        self._SolvePhysicsProblem()
        self._SetFunctionalWeights()
        self._SolveAdjointProblem()
        self.OutputSolutionStep()

    def _RunStageSolutionLoop(self):
        """
        This method executes a single physics solution of the Topology Optimization problem loop
        N.B.: must be called after the creation of the fluid-transport model part
        problem_stage = 1: Transport solution
        problem_stage = 2: Adjoint Transport solution
        """
        self._PrepareSettings()
        problem_stage = self.GetTopologyOptimizationStage()
        if ((problem_stage == 1)): 
            print("\n--|PHYSICS SOLUTION|")
        elif (problem_stage == 2):
            print("\n--|ADJOINT SOLUTION|")
        else:   
            print("--| UNKNOWN SOLUTION |")
        top_opt_stage_str = self.topology_optimization_stage_str
        print("--|" + top_opt_stage_str + "| START SOLUTION LOOP")
        while self.KeepAdvancingSolutionLoop():
            print("--|" + top_opt_stage_str + "| ADVANCE TIME")
            self.time = self._AdvanceTime()
            print("--|" + top_opt_stage_str + "| INITIALIZE SOLUTION STEP")
            self.InitializeSolutionStep()
            print("--|" + top_opt_stage_str + "| UPDATE PHYSICS PARAMETERS STEP")
            self.UpdatePhysicsParametersVariables()
            print("--|" + top_opt_stage_str + "| PREDICT")
            self._GetSolver().Predict()
            print("--|" + top_opt_stage_str + "| SOLVE SOLUTION STEP")
            is_converged = self._GetSolver().SolveSolutionStep()
            print("--|" + top_opt_stage_str + "| CHECK CONVERGENCE: skipped, it does not work! Why?")
            # self.__CheckIfSolveSolutionStepReturnsAValue(is_converged)
            print("--|" + top_opt_stage_str + "| FINALIZE SOLUTION STEP")
            self.FinalizeSolutionStep()
        print("--|" + top_opt_stage_str + "| END SOLUTION LOOP")

    def _SolvePhysicsProblem(self):
        self._SetTopologyOptimizationStage(1)
        self._RunStageSolutionLoop()

    def _SolveAdjointProblem(self):
        self._SetTopologyOptimizationStage(2)
        self._RunStageSolutionLoop()

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
    
    def _GetSolver(self, force_adjoint = False):
        """
        This method returns solver in use for the current topology optimization phase
        """
        if not hasattr(self, 'physics_solver'):
                self.physics_solver = self._CreateSolver()
        if not hasattr(self, 'adjoint_solver'):
                self.adjoint_solver = self._CreateSolver(True)
        self._solver = self._GetTopologyOptimizationStageSolver(force_adjoint)
        return self._solver
    
    def _GetPhysicsSolver(self):
        """
        This method returns the Transport Solver
        """
        if not hasattr(self, 'physics_solver'):
                self.physics_solver = self._CreateSolver(False)
        return self.physics_solver
    
    def _GetAdjointSolver(self):
        """
        This method returns the Adjoint Transport Solver
        """
        if not hasattr(self, 'adjoint_solver'):
                self.adjoint_solver = self._CreateSolver(True)
        return self.adjoint_solver
    
    def _GetTopologyOptimizationStageSolver(self, force_adjont = False):
        """
        This methods returns the current topology optimization stage solver
        iF force_adjoint --> return adjoint_solver
        If topology_optimization_stage != 2 (<=> EVERYTHING BUT NOT ADJ STAGE) --> return physics_solver
        If topology_optimization_stage == 2 (<=>  ADJ STAGE) --> return adjoint_solver
        """
        if (self.IsAdjointStage() or (force_adjont)): # ADJ
            return self.adjoint_solver
        else: # NS
            return self.physics_solver
        
    def IsInitializeStage(self):
        """
        This method returns TRUE the topology optimization stage is 0 <=> INIT
        """ 
        return self.CheckOptimizationStage(0)
    
    def IsPhysicsStage(self):
        """
        This method returns TRUE the topology optimization stage is 1 <=> NS
        """
        return self.CheckOptimizationStage(1)

    def IsAdjointStage(self):
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
        # self._SetFunctionalWeights()
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
        self.PreparePhysicsSolver()
        self.PrepareAdjointSolver()
    
    def PreparePhysicsSolver(self):
        """This method prepares the Navier-Stokes primal problem Solver in the AnalysisStage 
        Usage: It is designed to be called ONCE, BEFORE the execution of the solution-loop
        Prepare Solver : ImportModelPart -> PrepareModelPart -> AddDofs
        """
        self._GetPhysicsSolver().ImportModelPart()
        self._GetPhysicsSolver().PrepareModelPart()
        self._GetPhysicsSolver().AddDofs()
    
    def PrepareAdjointSolver(self):
        """This method prepares the Adjoint Navier-Stokes Solver in the AnalysisStage 
        Usage: It is designed to be called ONCE, BEFORE the execution of the solution-loop
        Prepare Solver : ImportModelPart -> PrepareModelPart -> AddDofs
        """
        # Modelers:
        self._GetAdjointSolver().ImportModelPart(self._GetPhysicsMainModelPartsList())
        self._GetAdjointSolver().PrepareModelPart()
        self._GetAdjointSolver().AddDofs()

    def InitializeSolvers(self):
        """This method initializes the NS and ADJ_NS Solvers in the AnalysisStage 
        Usage: It is designed to be called ONCE, BEFORE the execution of the solution-loop
        """
        self.InitializePhysicsSolver()
        self.InitializeAdjointSolver()

    def InitializePhysicsSolver(self):
        """This method initializes the NS Solver in the AnalysisStage 
        Usage: It is designed to be called ONCE, BEFORE the execution of the solution-loop
        """
        self._GetPhysicsSolver().Initialize()

    def InitializeAdjointSolver(self):
        """This method initializes the ADJ_NS Solver in the AnalysisStage 
        Usage: It is designed to be called ONCE, BEFORE the execution of the solution-loop
        """
        self._GetAdjointSolver().Initialize()

    def _SetFunctionalWeights(self, weights=[]):
        functional_weights = np.zeros(30)
        self._GetComputingModelPart().ProcessInfo.SetValue(KratosMultiphysics.FUNCTIONAL_WEIGHTS, functional_weights)
    
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
        if self.IsPhysicsStage(): # T
            return self.time < self.end_time
        elif self.IsAdjointStage(): #ADJ
            return self.time > self.start_time
        else:
            print("--|" + self.topology_optimization_stage_str + "| Time check outside T or ADJ_T solution")
            return False
        
    def _AdvanceTime(self):
        """
        This methods hadls the time for the Topology Optimization problems, currently its purpos it's just to make thinks work
        """
        time = super()._AdvanceTime()
        if self.IsAdjointStage():
            time = 0.0
        return time
    
    def _PrepareSettings(self):
        convention_diffusion_settings = self._GetComputingModelPart().ProcessInfo.GetValue(KratosMultiphysics.CONVECTION_DIFFUSION_SETTINGS)
        if (self.IsPhysicsStage()):
            convention_diffusion_settings.SetUnknownVariable(KratosMultiphysics.KratosGlobals.GetVariable("TEMPERATURE"))
            convention_diffusion_settings.SetVolumeSourceVariable(KratosMultiphysics.KratosGlobals.GetVariable("HEAT_FLUX"))
            convention_diffusion_settings.SetSurfaceSourceVariable(KratosMultiphysics.KratosGlobals.GetVariable("FACE_HEAT_FLUX"))
        elif (self.IsAdjointStage()):
            convention_diffusion_settings.SetUnknownVariable(KratosMultiphysics.KratosGlobals.GetVariable("TEMPERATURE_ADJ"))
            convention_diffusion_settings.SetVolumeSourceVariable(KratosMultiphysics.KratosGlobals.GetVariable("HEAT_FLUX_ADJ"))
            convention_diffusion_settings.SetSurfaceSourceVariable(KratosMultiphysics.KratosGlobals.GetVariable("FACE_HEAT_FLUX_ADJ"))
        else:
            print("\n ERROR: _PrepareSettings in the wrong topology optimization stage.\n")

    def GetTopologyOptimizationStage(self):
        """
        This method returns the topology optimization stage
        """
        return self.topology_optimization_stage

    def _GetPhysicsSolver(self):
        """
        This method returns the Navier-Stokes Solver
        """
        if not hasattr(self, 'physics_solver'):
                self.physics_solver = self._CreateSolver(False)
        return self.physics_solver
    
    def _GetAdjointSolver(self):
        """
        This method returns the Adjoint Navier-Stokes Solver
        """
        if not hasattr(self, 'adjoint_solver'):
                self.adjoint_solver = self._CreateSolver(True)
        return self.adjoint_solver

    def _GetPhysicsMainModelPartsList(self):
        return [self._GetPhysicsSolver().main_model_part]
    
    def _GetMainModelPart(self):
        """
        This method returns the main model part for the current solver
        """    
        return self._GetSolver().GetMainModelPart()
    
    def UpdatePhysicsParametersVariables(self):
        self._UpdateConductivityVariable()
        self._UpdateDecayVariable()
        self._UpdateConvectionCoefficientVariable()

    def _UpdateConductivityVariable(self):
        self._GetSolver()._UpdateConductivityVariable(self.conductivity)

    def _UpdateDecayVariable(self):
        self._GetSolver()._UpdateDecayVariable(self.decay)

    def _UpdateConvectionCoefficientVariable(self):
        self._GetSolver()._UpdateConvectionCoefficientVariable(self.convection_coefficient)

    def _InitializeTopologyOptimizationStepPhysicsSolution(self):
        """
        This method Re-Initializes the physical values in order to be consistent with the solution of the current optimization step
        """
        print("\n--|INITIALIZE OPTIMIZATION STEP PHYSICS SOLUTION|")
        self._SetTopologyOptimizationStage(0)
        if (not self.first_iteration):
                self._ReInitializePhysics()
        self._UpdatePhysicsParameters()

    def _ReInitializePhysics(self):
        """
        This method Re-Initializes the physical values in order to be consistent with the solution of the current optimization step
        """
        print("--|" + self.topology_optimization_stage_str + "| RE-INITIALIZE PHYSICS")
        self._GetComputingModelPart().ProcessInfo.SetValue(KratosMultiphysics.TIME, 0.0)

    def _UpdatePhysicsParameters(self):
        self._UpdateConductivity()
        self._UpdateDecay()
        self._UpdateConvectionCoefficient()

    def _UpdateConductivity(self):
        """
        This method handles the conductivity update.
        """
        print("--|" + self.topology_optimization_stage_str + "| UPDATE CONDUCTIVITY")
        self.conductivity = self._ComputeConductivity()
        self._UpdateConductivityVariable()

    def _UpdateDecay(self):
        """
        This method handles the decay update.
        """
        print("--|" + self.topology_optimization_stage_str + "| UPDATE DECAY")
        self.decay = self._ComputeDecay()
        self._UpdateDecayVariable()

    def _UpdateConvectionCoefficient(self):
        """
        This method handles the condcutivity update.
        """
        print("--|" + self.topology_optimization_stage_str + "| UPDATE CONVECTION COEFFICIENT")
        self.convection_coefficient = self._ComputeConvectionCoefficient()
        self._UpdateConvectionCoefficientVariable()
    
    def _ComputeConductivity(self):
        conductivity = self.conductivity_parameters[0]
        return conductivity
    
    def _ComputeDecay(self):
        decay = self.decay_parameters[0]
        return decay
    
    def _ComputeConvectionCoefficient(self):
        convection_coeff = self.convection_coefficient_parameters[0]
        return convection_coeff
    
    
    
    










