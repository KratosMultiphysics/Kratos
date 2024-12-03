from sys import argv

import nlopt #import the module for the optimizers implementation in python

import KratosMultiphysics as Kratos
# Import applications
import KratosMultiphysics.FluidDynamicsApplication as KratosCFD

from KratosMultiphysics.FluidDynamicsApplication.fluid_dynamics_analysis import FluidDynamicsAnalysis
from KratosMultiphysics.analysis_stage import AnalysisStage
from KratosMultiphysics.FluidDynamicsApplication import fluid_topology_optimization_solver

class FluidTopologyOptimizationAnalysis(FluidDynamicsAnalysis):
    def __init__(self,model,parameters):
        self.topology_optimization_stage = 0
        self.topology_optimization_stage_str = "INIT"
        super().__init__(model,parameters) 
        self._CreateTopologyOptimizationSolvers() # currently it is a useless method 
        self.max_it = 1
        

    def _CreateTopologyOptimizationSolvers(self):
        self.NS_solver = fluid_topology_optimization_solver.CreateSolver(self.model, self.project_parameters)
        self.ADJ_NS_solver = fluid_topology_optimization_solver.CreateSolver(self.model, self.project_parameters, True)

    def _CreateSolver(self, isAdjointSolver = False):
        return fluid_topology_optimization_solver.CreateSolver(self.model, self.project_parameters, isAdjointSolver)
    
    def _GetSolver(self, force_adjoint = False):
        if not hasattr(self, 'NS_solver'):
                self.NS_solver = self._CreateSolver()
        if not hasattr(self, 'ADJ_NS_solver'):
                self.ADJ_NS_solver = self._CreateSolver(True)
        self._solver = self._GetTopologyOptimizationStageSolver(force_adjoint)
        return self._solver
    
    def _GetNavierStokesSolver(self):
        if not hasattr(self, 'ADJ_NS_solver'):
                self.NS_solver = self._CreateSolver(False)
        return self.NS_solver
    
    def _GetAdjointNavierStokesSolver(self):
        if not hasattr(self, 'ADJ_NS_solver'):
                self.ADJ_NS_solver = self._CreateSolver(True)
        return self.ADJ_NS_solver
    
    def _GetTopologyOptimizationStageSolver(self, force_adjont = False):
        """
        This methods returns the current topology optimization stage solver
        iF force_adjoint --> return ADJ_NS_solver
        If topology_optimization_stage != 2 (<=> EVERYTHING BUT NOT ADJ STAGE) --> return NS_solver
        If topology_optimization_stage == 2 (<=>  ADJ STAGE) --> return ADJ_NS_solver
        """
        if (self.IsAdjointNavierStokesStage() or (force_adjont)): # ADJ
            return self.ADJ_NS_solver
        else: # NS
            return self.NS_solver

    def RunSolutionLoop(self):
        """
        This function executes the solution loop of the Fluid Topology Optimization Analysis
        """
        it = 0
        while (it < self.max_it):
            ## Advance Iteration
            it = it+1
            print("\n--------------------------------------------------------")
            print(  "--| FLUID TOPOLOGY OPTIMIZATION SOLUTION LOOP. IT:", it)
            print("--------------------------------------------------------")
            self._UpdateResistance()
            ## NAVIER-STOKES SOLUTION
            self._RunStageSolutionLoop(1) # NAVIER-STOKES PROBLEM SOLUTION
            ## ADJOINT NAVIER-STOKES SOLUTION
            self._RunStageSolutionLoop(2) # ADJOINT NAVIER-STOKES PROBLEM SOLUTION
            ## DOMAIN DESIGN UPDATE
            self._UpdateDomainDesign()  # DESIGN PARAMETER UPDATE
            print("--------------------------------------------------------\n")

    def _RunStageSolutionLoop(self, problem_stage):
        """
        This function executes a single physics solution loop of the Topology Optimization problem.
        N.B.: must be called after the creation of the fluid model part
        problem_stage = 1: Navier Stokes solution
        problem_stage = 2: Adjoint Navier-Stokes solution
        """
        if (problem_stage == 1): # NS
            print("\n--|NAVIER-STOKES SOLUTION|")
        elif (problem_stage == 2): # ADJ
            print("\n--|ADJOINT NAVIER-STOKES SOLUTION|")
        else:
            print("--| UNKNOWN SOLUTION |")
        self._SetTopologyOptimizationStage(problem_stage)
        top_opt_stage_str = self.topology_optimization_stage_str
        print("--|" + top_opt_stage_str + "| START SOLUTION LOOP")
        while self.KeepAdvancingSolutionLoop():
            print("--|" + top_opt_stage_str + "| ADVANCE TIME")
            self.time = self._AdvanceTime()
            print("--|" + top_opt_stage_str + "| INITIALIZE SOLUTION STEP")
            self.InitializeSolutionStep()
            print("--|" + top_opt_stage_str + "| PREDICT")
            self._GetSolver().Predict()
            print("--|" + top_opt_stage_str + "| SOLVE SOLUTION STEP")
            is_converged = self._GetSolver().SolveSolutionStep()
            print("--|" + top_opt_stage_str + "| CHECK CONVERGENCE: skipped, it does not work! Why?")
            # self.__CheckIfSolveSolutionStepReturnsAValue(is_converged)
            print("--|" + top_opt_stage_str + "| FINALIZE SOLUTION STEP")
            self.FinalizeSolutionStep()
            print("--|" + top_opt_stage_str + "| PRINT SOLUTION STEP OUTPUT")
            self.OutputSolutionStep()
        print("--|" + top_opt_stage_str + "| END SOLUTION LOOP")

    def _SetTopologyOptimizationStage(self, problem_stage):
        """
        This method set the topology optimization problem phase:
        0: initialization, 1: NS solution, 2: ADJ NS solution, 3: Optimization Update
        """
        self.topology_optimization_stage = problem_stage 
        if self.CheckOptimizationStage(0):
            self.topology_optimization_stage_str = "INIT"
        elif self.CheckOptimizationStage(1):
            self.topology_optimization_stage_str = "NS"
        elif self.CheckOptimizationStage(2):
            self.topology_optimization_stage_str = "ADJ"
        elif self.CheckOptimizationStage(3):
            self.topology_optimization_stage_str = "OPT"
        else:
            self.topology_optimization_stage_str = "ERROR"
        self._GetSolver().GetComputingModelPart().ProcessInfo.SetValue(KratosCFD.FLUID_TOP_OPT_PROBLEM_STAGE, self.topology_optimization_stage)

    def Initialize(self):
        """This function initializes the FluidTopologyOptimizationAnalysis
        Usage: It is designed to be called ONCE, BEFORE the execution of the solution-loop
        """
        # Modelers:
        self._CreateModelers()
        self._ModelersSetupGeometryModel()
        self._ModelersPrepareGeometryModel()
        self._ModelersSetupModelPart()

        # Prepare Solvers : ImportModelPart -> PrepareModelPart -> AddDofs)
        self.PrepareSolvers()

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
        if self._GetSolver().GetComputingModelPart().ProcessInfo[Kratos.IS_RESTARTED]:
            self.time = self._GetSolver().GetComputingModelPart().ProcessInfo[Kratos.TIME]
        else:
            self.time = self.project_parameters["problem_data"]["start_time"].GetDouble()
            self._GetSolver().GetComputingModelPart().ProcessInfo[Kratos.TIME] = self.time

        self.start_time = self.time
        self.end_time = self.project_parameters["problem_data"]["end_time"].GetDouble()

        ## If the echo level is high enough, print the complete list of settings used to run the simulation
        if self.echo_level > 1:
            with open("ProjectParametersOutput.json", 'w') as parameter_output_file:
                parameter_output_file.write(self.project_parameters.PrettyPrintJsonString())

        # Print Start Analysis
        Kratos.Logger.PrintInfo(self._GetSimulationName(), "Analysis -START- ")

    def PrepareSolvers(self):
        """This function prepares the NS and ADJ_NS Solvers in the AnalysisStage 
        Usage: It is designed to be called ONCE, BEFORE the execution of the solution-loop
        Prepare Solvers : ImportModelPart -> PrepareModelPart -> AddDofs
        """
        self.PrepareNavierStokesSolver()
        self.PrepareAdjointNavierStokesSolver()
    
    def PrepareNavierStokesSolver(self):
        """This function prepares the Navier-Stokes primal problem Solver in the AnalysisStage 
        Usage: It is designed to be called ONCE, BEFORE the execution of the solution-loop
        Prepare Solver : ImportModelPart -> PrepareModelPart -> AddDofs
        """
        self._GetNavierStokesSolver().ImportModelPart()
        self._GetNavierStokesSolver().PrepareModelPart()
        self._GetNavierStokesSolver().AddDofs()
    
    def PrepareAdjointNavierStokesSolver(self):
        """This function prepares the Adjoint Navier-Stokes Solver in the AnalysisStage 
        Usage: It is designed to be called ONCE, BEFORE the execution of the solution-loop
        Prepare Solver : ImportModelPart -> PrepareModelPart -> AddDofs
        """
        # Modelers:
        self._GetAdjointNavierStokesSolver().ImportModelPart()
        self._GetAdjointNavierStokesSolver().PrepareModelPart()
        self._GetAdjointNavierStokesSolver().AddDofs()

    def InitializeSolvers(self):
        """This function initializes the NS and ADJ_NS Solvers in the AnalysisStage 
        Usage: It is designed to be called ONCE, BEFORE the execution of the solution-loop
        """
        self.InitializeNavierStokesSolver()
        self.InitializeAdjointNavierStokesSolver()

    def InitializeNavierStokesSolver(self):
        """This function initializes the NS Solver in the AnalysisStage 
        Usage: It is designed to be called ONCE, BEFORE the execution of the solution-loop
        """
        self._GetNavierStokesSolver().Initialize()

    def InitializeAdjointNavierStokesSolver(self):
        """This function initializes the ADJ_NS Solver in the AnalysisStage 
        Usage: It is designed to be called ONCE, BEFORE the execution of the solution-loop
        """
        self._GetAdjointNavierStokesSolver().Initialize()

    def Check(self):
        """This function checks the AnalysisStage
        Keyword arguments:
        self -- It signifies an instance of a class.
        """
        # Checking solver
        self._GetNavierStokesSolver().Check()
        self._GetAdjointNavierStokesSolver().Check()

        # Checking processes
        for process in self._GetListOfProcesses():
            process.Check()

    def __CreateListOfProcesses(self):
        """This function creates the processes and the output-processes
        """
        order_processes_initialization = self._GetOrderOfProcessesInitialization()
        self._list_of_processes        = self._CreateProcesses("processes", order_processes_initialization)
        deprecated_output_processes    = self._CheckDeprecatedOutputProcesses(self._list_of_processes)
        order_processes_initialization = self._GetOrderOfOutputProcessesInitialization()
        self._list_of_output_processes = self._CreateProcesses("output_processes", order_processes_initialization)
        self._list_of_processes.extend(self._list_of_output_processes) # Adding the output processes to the regular processes
        self._list_of_output_processes.extend(deprecated_output_processes)

    def _UpdateDomainDesign(self):
        """
        This method will handle the design parameter update.
        1) functional derivatives evaluation
        2) optimizer (MMA)?
        """
        self._SetTopologyOptimizationStage(3)
        print("--|" + self.topology_optimization_stage_str + "| DESIGN PARAMETER UPDATE NOT YET IMPLEMENTED!")
        self._EvaluateFunctionalDerivatives()
        self._SolveOptimizer()  

    def _EvaluateFunctionalValueAndDerivatives(self):
        """
        This method is used to evaluate the functional value and derivatives w.r.t. the design parameter
        """
        self._EvaluateFunctionalValue()
        self._EvaluateFunctionalDerivatives()
        dumb = 0
        return dumb
    
    def _EvaluateFunctionalValue(self):
        """
        This method is used to evaluate the functional value
        """
        dumb = 0
        return dumb
    
    def _EvaluateFunctionalDerivatives(self):
        """
        This method is used to evaluate the functional derivatives w.r.t the design parameter
        """
        dumb = 0
        return dumb
    
    def _SolveOptimizer(self):
        """
        This method is used to update the design parameter using an optimizer
        """
        dumb = 0
        return dumb
    
    def _UpdateResistance(self):
        """
        This method will handle the resistance update, that will be based on the value of the design parameter.
        Now it is used to "play" with the resistance
        """
        dumb = 0
        return dumb

    
    def KeepAdvancingSolutionLoop(self):
        """This function specifies the stopping criteria for breaking the solution loop
        It can be overridden by derived classes
        """
        if self.IsNavierStokesStage(): # NS
            return self.time < self.end_time
        elif self.IsAdjointNavierStokesStage(): #ADJ
            return self.time > self.start_time
            # print("--|" + self.topology_optimization_stage_str + "| AVOIDED SOLUTION")
            # return False
        else:
            print("--|" + self.topology_optimization_stage_str + "| Time check outside NS or ADJ_NS solution")
            return False
        
    def PrintAnalysisStageProgressInformation(self):
        Kratos.Logger.PrintInfo(self._GetSimulationName(), "TOTAL STEP: ", self._GetSolver().GetComputingModelPart().ProcessInfo[Kratos.STEP])
        if self.IsNavierStokesStage(): # NS
            Kratos.Logger.PrintInfo(self._GetSimulationName(), "NS STEP: ", self._GetSolver().GetComputingModelPart().ProcessInfo[KratosCFD.FLUID_TOP_OPT_NS_STEP])
        elif self.IsAdjointNavierStokesStage(): # ADJ
            Kratos.Logger.PrintInfo(self._GetSimulationName(), "ADJ NS STEP: ", self._GetSolver().GetComputingModelPart().ProcessInfo[KratosCFD.FLUID_TOP_OPT_ADJ_NS_STEP])
        else:
            Kratos.Logger.PrintInfo(self._GetSimulationName(), "Invalid Value of the Topology Optimization Stage. TOP_OPT_STAGE: ", self.topology_optimization_stage)
        Kratos.Logger.PrintInfo(self._GetSimulationName(), "TIME: ", self.time)

    def _GetSimulationName(self):
        return "--|" + self.topology_optimization_stage_str + "| Solution Analysis"
    
    def _GetOrderOfProcessesInitialization(self):
        # The list of processes will contain a list with each individual process already constructed (boundary conditions, initial conditions and gravity)
        # Note 1: gravity is constructed first. Outlet process might need its information.
        # Note 2: initial conditions are constructed before BCs. Otherwise, they may overwrite the BCs information.
        # Note 3: ADJ NS boundary conditions are constructed after the NS boundary condition, in they could require them.
        return ["gravity",
                "initial_conditions_process_list",
                "boundary_conditions_process_list",
                "adjoint_boundary_conditions_process_list"
                "auxiliar_process_list"]

    def Run(self):
        """
        This function executes the entire Fluid Topology Optimization Stage
        """
        self.Initialize()
        self.RunSolutionLoop()
        self.Finalize()

    def ApplyBoundaryConditions(self):
        """here the boundary conditions is applied, by calling the InitializeSolutionStep function of the processes"""

        for process in self._GetListOfProcesses():
            process.ExecuteInitializeSolutionStep()

        # bc_model_part = ...
        # for node in bc_model_part.Nodes:
            #if(node.IsFixed(Kratos.VELOCITY_X) and )

        #other operations as needed

    def GetOptimizationStage(self):
        return self.topology_optimization_stage

    def CheckOptimizationStage(self, check):
        return (self.GetOptimizationStage() == check)
    
    def IsInitializeStage(self):
        return self.CheckOptimizationStage(0)
    
    def IsNavierStokesStage(self):
        return self.CheckOptimizationStage(1)

    def IsAdjointNavierStokesStage(self):
        return self.CheckOptimizationStage(2)
    
    def IsOptimizerStage(self):
        return self.CheckOptimizationStage(3)
    
    def _AdvanceTime(self):
        time = super()._AdvanceTime()
        if self.IsAdjointNavierStokesStage():
            time = 0.0
        return time







        