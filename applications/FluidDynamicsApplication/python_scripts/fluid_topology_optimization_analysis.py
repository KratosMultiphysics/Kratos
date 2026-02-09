## IMPORT
# Import Libraries
from sys import argv
import numpy as np #import the numpy library
import scipy as sp #import the scipy library
from scipy.spatial import KDTree
from scipy.sparse import dok_matrix, lil_matrix
import mmapy as MMA # import the MMA subroutines python library: https://github.com/arjendeetman/GCMMA-MMA-Python
import time
import shutil
from pathlib import Path

# Import Kratos
import KratosMultiphysics as KratosMultiphysics
from KratosMultiphysics.kratos_utilities import IssueDeprecationWarning

# MPI utilities
# Auxiliary function to check the parallel type at runtime
#TODO: Delete this once we come up with the final factory-based design
def _CheckIsDistributed():
    if KratosMultiphysics.ParallelEnvironment.HasDataCommunicator("World"):
        world_data_comm = KratosMultiphysics.ParallelEnvironment.GetDataCommunicator("World")
        return world_data_comm.IsDistributed()
    else:
        return False
# If required, import parallel applications and modules
if _CheckIsDistributed():
    import KratosMultiphysics.mpi as KratosMPI
    import KratosMultiphysics.MetisApplication as KratosMetis
    import KratosMultiphysics.TrilinosApplication as KratosTrilinos
    import KratosMultiphysics.mpi.distributed_import_model_part_utility as distributed_import_model_part_utility
# Importing factories
if _CheckIsDistributed():
    import KratosMultiphysics.TrilinosApplication.trilinos_linear_solver_factory as linear_solver_factory
else:
    import KratosMultiphysics.python_linear_solver_factory as linear_solver_factory
    import KratosMultiphysics.base_convergence_criteria_factory as convergence_criteria_factory

from KratosMultiphysics import DataCommunicator

# Import Kratos Applications
import KratosMultiphysics.FluidDynamicsApplication as KratosCFD
import KratosMultiphysics.ConvectionDiffusionApplication as KratosCD
import KratosMultiphysics.MeshingApplication as KratosMMG
# Import Kratos Analysis and Solvers
from KratosMultiphysics.FluidDynamicsApplication.fluid_dynamics_analysis import FluidDynamicsAnalysis
from KratosMultiphysics.FluidDynamicsApplication import trilinos_fluid_topology_optimization_solver
from KratosMultiphysics.FluidDynamicsApplication import fluid_topology_optimization_solver
from KratosMultiphysics.ConvectionDiffusionApplication.apply_topology_optimization_pde_filter_process import ApplyTopologyOptimizationPdeFilterProcess
# Import Kratos Processes
from KratosMultiphysics import ComputeNodalGradientProcess

class FluidTopologyOptimizationAnalysis(FluidDynamicsAnalysis):
    def _ReadOptimizationParameters(self):
        self.project_parameters_adjoint = self.project_parameters.Clone()
        if (self.project_parameters.Has("optimization_parameters_file_name")):
            self.optimization_parameters_file = self.project_parameters["optimization_parameters_file_name"].GetString()
        else:
            self.optimization_parameters_file = "OptimizationParameters.json"
        with open(self.optimization_parameters_file, 'r') as parameter_file:
            self.optimization_parameters = KratosMultiphysics.Parameters(parameter_file.read())
        self.physics_parameters_settings = self.optimization_parameters["physics_parameters_settings"]
        self.optimization_settings       = self.optimization_parameters["optimization_settings"]
        self.ValidateOptimizationParameters()

    def _SetTopologyOptimizationName(self):
        self.topology_optimization_name = "FLUID"

    def _CreateTopologyOptimizationSolvers(self):
        """
        This method creates the NS and ADJ_NS solvers
        """
        if self.IsMpiParallelism():
            self.physics_solver = trilinos_fluid_topology_optimization_solver.CreateSolver(self.model, self.project_parameters)
            self.adjoint_solver = trilinos_fluid_topology_optimization_solver.CreateSolver(self.model, self.project_parameters_adjoint, isAdjointSolver=True)
        else:
            self.physics_solver = fluid_topology_optimization_solver.CreateSolver(self.model, self.project_parameters)
            self.adjoint_solver = fluid_topology_optimization_solver.CreateSolver(self.model, self.project_parameters_adjoint, isAdjointSolver=True)

    def _CreateSolver(self, isAdjointSolver = False):
        """
        This method creates a solver
        isAdjointSolver == False --> physics_solver
        isAdjointSolver == True  --> adjoint_solver
        """
        if self.IsMpiParallelism():
            if (not isAdjointSolver):
                return trilinos_fluid_topology_optimization_solver.CreateSolver(self.model, self.project_parameters)
            else:
                return trilinos_fluid_topology_optimization_solver.CreateSolver(self.model, self.project_parameters_adjoint, isAdjointSolver=True)
        else:
            if (not isAdjointSolver):
                return fluid_topology_optimization_solver.CreateSolver(self.model, self.project_parameters, isAdjointSolver)
            else:
                return fluid_topology_optimization_solver.CreateSolver(self.model, self.project_parameters_adjoint, isAdjointSolver=True)
    
    def _GetSolver(self, force_adjoint = False):
        """
        This method returns solver in use for the current topology optimization phase
        """
        if not hasattr(self, 'physics_solver'):
                self.physics_solver = self._CreateSolver()
        if not hasattr(self, 'adjoint_solver'):
                self.adjoint_solver = self._CreateSolver(isAdjointSolver=True)
        self._solver = self._GetTopologyOptimizationStageSolver(force_adjoint)
        return self._solver
    
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
    
    def _GetTopologyOptimizationStageSolver(self, force_adjoint = False):
        """
        This methods returns the current topology optimization stage solver
        iF force_adjoint --> return adjoint_solver
        If topology_optimization_stage != 2 (<=> EVERYTHING BUT NOT ADJ STAGE) --> return physics_solver
        If topology_optimization_stage == 2 (<=>  ADJ STAGE) --> return adjoint_solver
        """
        if (self.IsAdjointStage() or (force_adjoint)): # ADJ
            return self.adjoint_solver
        else: # NS
            return self.physics_solver
        
    def _SetMinMaxIt(self):
        """
        This method sets the minimum & maximum iteration for the topology optimization process
        """
        iterations_settings = self.optimization_settings["min_max_iterations"]
        self.min_it = iterations_settings[0].GetInt()
        self.max_it = iterations_settings[1].GetInt()
        if (self.min_it > self.max_it):
            self.MpiPrint("\n!!!WARNING: wrong initialization of the min & max number of iterations\n")
        
    def _CreateProcesses(self, parameter_name, initialization_order):
        """Create a list of Processes
        This method is TEMPORARY to not break existing code
        It will be removed in the future
        """
        processes_params = self.project_parameters[parameter_name]
        solution_output_path = self.optimization_parameters["optimization_settings"]["solution_output_settings"]["output_path"].GetString()
        if parameter_name == "optimization_output_processes":
            for name, value in processes_params.items():
                self._SetOutputProcessPath(value, add_path=solution_output_path+"/optimization")      
        elif parameter_name == "time_output_processes":
            for name, value in processes_params.items():
                self._SetOutputProcessPath(value, add_path=solution_output_path+"/time")
            
        list_of_processes = super(FluidDynamicsAnalysis, self)._CreateProcesses(parameter_name, initialization_order)

        if parameter_name == "processes":
            processes_block_names = ["gravity", "initial_conditions_process_list", "boundary_conditions_process_list", "auxiliar_process_list"]
            if len(list_of_processes) == 0: # Processes are given in the old format (or no processes are specified)
                for process_name in processes_block_names:
                    if self.project_parameters.Has(process_name):
                        info_msg  = "Using the old way to create the processes, this was removed!\n"
                        info_msg += "Refer to \"https://github.com/KratosMultiphysics/Kratos/wiki/Common-"
                        info_msg += "Python-Interface-of-Applications-for-Users#analysisstage-usage\" "
                        info_msg += "for a description of the new format"
                        raise Exception("FluidTransportTopologyOptimizationAnalysis: " + info_msg)
            else: # Processes are given in the new format
                for process_name in processes_block_names:
                    if (self.project_parameters.Has(process_name) is True):
                        raise Exception("Mixing of process initialization is not alowed!")
        elif parameter_name == "optimization_output_processes":
            if self.project_parameters.Has("output_configuration"):
                info_msg  = "Using the old way to create the gid-output, this was removed!\n"
                info_msg += "Refer to \"https://github.com/KratosMultiphysics/Kratos/wiki/Common-"
                info_msg += "Python-Interface-of-Applications-for-Users#analysisstage-usage\" "
                info_msg += "for a description of the new format"
                raise Exception("FluidTransportTopologyOptimizationAnalysis: " + info_msg)
        elif parameter_name == "time_output_processes":
            if self.project_parameters.Has("output_configuration"):
                info_msg  = "Using the old way to create the gid-output, this was removed!\n"
                info_msg += "Refer to \"https://github.com/KratosMultiphysics/Kratos/wiki/Common-"
                info_msg += "Python-Interface-of-Applications-for-Users#analysisstage-usage\" "
                info_msg += "for a description of the new format"
                raise Exception("FluidTransportTopologyOptimizationAnalysis: " + info_msg)
        elif parameter_name == "output_processes":
            info_msg  =  "Using the standard analysis way of setting outputs, not accepatable for Fluid-Transport topology optimization analysis\n"
            info_msg  += "Replace it with \"optimization_output_processes\"\n"
            info_msg  += "If time dependent output is needed add a list named \"time_output_processes\"\n"
            raise Exception("FluidTransportTopologyOptimizationAnalysis: " + info_msg)
        else:
            raise NameError("wrong parameter name")

        return list_of_processes

    def __CreateListOfProcesses(self):
        """This method creates the processes and the output-processes
        """
        order_processes_initialization = self._GetOrderOfProcessesInitialization()
        self._list_of_processes        = self._CreateProcesses("processes", order_processes_initialization)
        deprecated_output_processes    = self._CheckDeprecatedOutputProcesses(self._list_of_processes)
        order_processes_initialization = self._GetOrderOfOutputProcessesInitialization()
        self._list_of_optimization_output_processes = self._CreateProcesses("optimization_output_processes", [])
        self._list_of_time_output_processes         = self._CreateProcesses("time_output_processes", [])
        self._list_of_output_processes = self._list_of_optimization_output_processes + self._list_of_time_output_processes
        self._list_of_processes.extend(self._list_of_output_processes) # Adding the output processes to the regular processes
        self._list_of_output_processes.extend(deprecated_output_processes)

    def _GetOrderOfProcessesInitialization(self):
        # The list of processes will contain a list with each individual process already constructed (boundary conditions, initial conditions and gravity)
        # Note 1: gravity is constructed first. Outlet process might need its information.
        # Note 2: initial conditions are constructed before BCs. Otherwise, they may overwrite the BCs information.
        # Note 3: ADJ NS boundary conditions are constructed after the NS boundary condition, in they could require them.
        return ["gravity",
                "initial_conditions_process_list",
                "boundary_conditions_process_list",
                "adjoint_boundary_conditions_process_list",
                "auxiliar_process_list"]

    def Run(self):
        """
        This method executes the entire TOPOLOGY Optimization Stage
        """
        self.Initialize()
        self.RunSolutionLoop()
        self.Finalize()

    def Initialize(self):
        """This method initializes the FluidTransportTopologyOptimizationAnalysis
        Usage: It is designed to be called ONCE, BEFORE the execution of the solution-loop
        """
        # Modelers:
        self._CreateModelers()
        self._ModelersSetupGeometryModel()
        self._ModelersPrepareGeometryModel()
        self._ModelersSetupModelPart()
        self.InitializeAnalysisTimeSettings()
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

        for process in self._GetListOfProcesses():
            process.ExecuteBeforeSolutionLoop()

        if self._GetComputingModelPart().ProcessInfo[KratosMultiphysics.IS_RESTARTED]:
            self.time = self._GetComputingModelPart().ProcessInfo[KratosMultiphysics.TIME]
        else:
            self.time = self.project_parameters["problem_data"]["start_time"].GetDouble()
            self._GetComputingModelPart().ProcessInfo[KratosMultiphysics.TIME] = self.time

        # Set Topology Optimization Stage: Initialize
        self._SetTopologyOptimizationStage(0)
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
        if self.IsMpiParallelism():
            self._GetAdjointSolver().ImportModelPart(model_parts=self._GetPhysicsMainModelPartsList(), physics_solver_distributed_model_part_importer=self._GetPhysicsSolverDistributedModelPartImporter())
        else:
            self._GetAdjointSolver().ImportModelPart(self._GetPhysicsMainModelPartsList())
        self._GetAdjointSolver().PrepareModelPart()
        self._GetAdjointSolver().AddDofs()

    def _SetOutputProcessPath(self, output_processes_settings, add_path=""):
        if output_processes_settings[0]["Parameters"].Has("output_path" ):
            base_output_path = output_processes_settings[0]["Parameters"]["output_path"].GetString()
        else:
            info_msg = "Defining an output process without settig the output path\n"
            raise Exception("FluidTransportTopologyOptimizationAnalysis: " + info_msg)
        output_processes_settings[0]["Parameters"]["output_path"].SetString(add_path + "/" + base_output_path)

    def _SetFunctionalWeights(self):
        # set future transport functionals to zero
        self._InitializeFunctionalWeights()
        self._InitializeFunctionalWeightsInTime()
        self._PrintFunctionalWeights()
        self.EvaluateTotalFunctional()

    def _InitializeFunctionalWeightsInTime(self):
        # the functional weights in time are used to solve the adjopint pde consistently, then in the functional evaluation the standard functional weigths (self.functional_weigths) are applied, since the time integration in performed autonomously.
        self.functional_weights_in_time = np.zeros((self.n_time_steps, self.n_functionals))
        self.InitializePhysicsFunctionalsWeightsInTime()

    def InitializePhysicsFunctionalsWeightsInTime(self):
        self._InitializeFluidFunctionalsWeightsInTime()

    def _InitializeFluidFunctionalsWeightsInTime(self):
        self._InitializeResistanceFunctionalWeightsInTime()
        self._InitializeStrainRateFunctionalWeightsInTime()
        self._InitializeVorticityFunctionalWeightsInTime()

    def _InitializeResistanceFunctionalWeightsInTime(self):
        resistance_functional_id = 0
        self._SetFunctionalWeightsInTime(resistance_functional_id, self.resistance_functional_time_info)

    def _InitializeStrainRateFunctionalWeightsInTime(self):
        strain_rate_functional_id = 1
        self._SetFunctionalWeightsInTime(strain_rate_functional_id, self.strain_rate_functional_time_info)

    def _InitializeVorticityFunctionalWeightsInTime(self):
        vorticity_functional_id = 2
        self._SetFunctionalWeightsInTime(vorticity_functional_id, self.vorticity_functional_time_info)

    def _SetFunctionalWeightsInTime(self, functional_id, time_info):
        """ Set the functional time-weights to 1.0 over the activation interval [start_time, end_time].
        start_time and end_time are continuous and each lies between two discrete time steps
        (the corresponding span is given by *_span_time_steps_ids and *_span_times).
        We map each boundary to one of the two steps by comparing it with mean(span_times):
          - start: pick the first id if start_time <= mean(start_span_times), else the second.
          - end  : pick the first id if end_time   <= mean(end_span_times),   else the second.
        Then we activate all time steps in [start_weights_id, end_weights_id] (inclusive).
        If start_weights_id < 0, clamp it to 0 (weights start from the first time step)."""
        start_time = time_info["times"][0]
        end_time   = time_info["times"][1]
        if (start_time <= np.mean(time_info["start_span_times"])):
            start_weights_id = time_info["start_span_time_steps_ids"][0]
        else:
            start_weights_id = time_info["start_span_time_steps_ids"][1]
        if (end_time <= np.mean(time_info["end_span_times"])):
            end_weights_id = time_info["end_span_time_steps_ids"][0]
        else:
            end_weights_id = time_info["end_span_time_steps_ids"][1]
        if start_weights_id < 0:
            start_weights_id = 0
        for time_step_id in range(start_weights_id, end_weights_id+1):
            self.functional_weights_in_time[time_step_id, functional_id] = 1.0

    def _PrintFunctionalWeights(self):
        self._PrintFunctionalWeightsPhysicsInfo()
        self.MpiPrint("--|" + self.topology_optimization_stage_str + "| REAL FUNCTIONAL WEIGHTS: " + str(self.functional_weights))
        self.MpiPrint(self.optimization_settings["optimization_problem_settings"]["functional_weights"].PrettyPrintJsonString())

    def _PrintFunctionalWeightsPhysicsInfo(self):
        self.MpiPrint("--|" + self.topology_optimization_stage_str + "| INITIAL FUNCTIONAL: " + str(self.initial_fluid_functional))
        self.MpiPrint("--|" + self.topology_optimization_stage_str + "| NORMALIZED FUNCTIONAL WEIGHTS: " + str(self.normalized_fluid_functional_weights))
    
    def _InitializeFunctionalWeights(self):
        "This method assumes that the '_ImportFunctionalWeights()' has already been called"
        if (self.functional_weights_imported):
            # get number of functionals
            self.n_fluid_functionals = len(self.normalized_fluid_functional_weights)
            self.n_functionals = self.n_fluid_functionals + 7 # 7 = transportfunctionals
            # initialize initial functionals vector container
            self.initial_fluid_functionals_values = np.zeros(self.n_fluid_functionals)
            self.initial_functionals_values = np.zeros(self.n_functionals)
            # initialize functionals vector container
            self.fluid_functionals = np.zeros(self.n_fluid_functionals)
            self.functionals = np.zeros(self.n_functionals)
            self.EvaluateFunctionals(print_functional=False)
            self._SetInitialFunctionals()
            self._SetNormalizationFunctionals()
            self.functional_weights = self._RescaleFunctionalWeightsByNormalizationValues()
        else:
            info_msg = "Calling '_InitializeFunctionalWeights' method before '_ImportFunctionalWeights()'"
            raise RuntimeError("FluidTopologyOptimizationAnalysis: " + info_msg)

    def _ImportFunctionalWeights(self):
        self.ImportPhysicsFunctionalWeights()
        self.functional_weights_imported = True

    def ImportPhysicsFunctionalWeights(self):
        self._ImportFluidFunctionalWeights()

    def _ImportFluidFunctionalWeights(self):
        fluid_functional_weights = [0.0, 0.0, 0.0]
        functional_weights_parameters = self.optimization_settings["optimization_problem_settings"]["functional_weights"]["fluid_functionals"]
        self.fluid_functional_normalization_strategy = functional_weights_parameters["normalization"]["type"].GetString()
        if self.fluid_functional_normalization_strategy == "custom":
            self.fluid_functional_normalization_value = functional_weights_parameters["normalization"]["value"].GetDouble()
        else:
            self.fluid_functional_normalization_strategy == "initial"
        fluid_functional_weights[0] = functional_weights_parameters["resistance"]["weight"].GetDouble()
        fluid_functional_weights[1] = functional_weights_parameters["strain_rate"]["weight"].GetDouble()
        fluid_functional_weights[2] = functional_weights_parameters["vorticity"]["weight"].GetDouble()
        # normalize weights
        self.normalized_fluid_functional_weights = self._NormalizeFunctionalWeights(np.asarray(fluid_functional_weights))
    
    def _SetInitialFunctionals(self):
        self.initial_fluid_functional = np.dot(self.normalized_fluid_functional_weights, self.initial_fluid_functionals_values)
        if (abs(self.initial_fluid_functional) < 1e-10):
            self.MpiPrint("[WARNING] Initial fluid functional is zero")

    def _SetNormalizationFunctionals(self):
        if self.fluid_functional_normalization_strategy == "initial":
            self.fluid_functional_normalization_value = self.initial_fluid_functional
        else: # custom value, already defined
            pass
        
    def _RescaleFunctionalWeightsByNormalizationValues(self):
        if (np.sum(np.abs(self.normalized_fluid_functional_weights)) < 1e-10):
            fluid_functional_weights = np.zeros(self.normalized_fluid_functional_weights.size)
        else:
            fluid_functional_weights  = self.normalized_fluid_functional_weights.copy()
            if (abs(self.initial_fluid_functional) > 1e-15):
                fluid_functional_weights /= abs(self.fluid_functional_normalization_value)
        return np.concatenate((fluid_functional_weights, np.zeros(self.n_functionals-self.n_fluid_functionals)))

    def _NormalizeFunctionalWeights(self, weights):
        weights_sum = np.sum(np.abs(weights))
        if (weights_sum > 1e-10):
            return np.asarray(weights)/weights_sum
        else:
            return np.asarray(weights)
    
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

    def InitializeAnalysisTimeSettings(self):
        self._InitializeAnalysisTimeCoefficient()
        self._InitializeAnalysisTimeStepsSettings()

    def _InitializeAnalysisTimeCoefficient(self):
        self._GetPhysicsSolver()._SetAnalysisTimeCoefficient()

    def _InitializeAnalysisTimeStepsSettings(self):    
        # Stepping and time settings
        if self.IsUnsteadySolution():
            self.start_time = self.project_parameters["problem_data"]["start_time"].GetDouble()
            self.end_time = self.project_parameters["problem_data"]["end_time"].GetDouble()
            self.time_step_base = self._GetSolver().settings["time_stepping"]["time_step"].GetDouble()
            self.n_time_steps = max(int(np.ceil((self.end_time-self.start_time) / self.time_step_base)), 1)
            self.simulation_times = np.zeros(self.n_time_steps+1)
            self.time_steps = np.zeros(self.n_time_steps)
            self.time_steps_integration_weights = np.zeros(self.n_time_steps) # weight to each time step to perform fast time integration of functionals
            for istep in range(self.n_time_steps):
                self.time_steps[istep] = self.time_step_base # assume equal time steps
                self.simulation_times[istep] = self.start_time + istep*self.time_step_base
            self.simulation_times[self.n_time_steps] = self.end_time
            # now build th eintegration weigths
            # the idea is that the ith time solution will be evaluated (i-1)th time step
            # functional at time_id=n multiplies the time_step_id=n-1
            # particular case: time_id=0, is not evaluated since it is not coupled with any time step. It is supposed that vel(t=t0)=0
            for istep in range(self.n_time_steps-1):
                self.time_steps_integration_weights[istep] = (self.time_steps[istep] + self.time_steps[istep+1]) / 2.0
            if (self.n_time_steps == 1):
                self.time_steps_integration_weights[self.n_time_steps-1] = self.time_steps[self.n_time_steps-1]
            else:
                self.time_steps_integration_weights[self.n_time_steps-1] = self.time_steps[self.n_time_steps-1] / 2.0
        else:
            self.start_time = self.project_parameters["problem_data"]["start_time"].GetDouble()
            self.end_time = self.start_time + 1e10
            self.time_step_base = 1.0000000001*1e10
            self.n_time_steps = 1
            self.time_steps = np.zeros(self.n_time_steps)
            self.time_steps_integration_weights = np.zeros(self.n_time_steps) # weight to each time step to perform fast time integration of functionals
            self.time_steps[0] = self.time_step_base
            self.time_steps_integration_weights[0] = 1.0
            

    def Check(self):
        """This method checks the AnalysisStage
        Keyword arguments:
        self -- It signifies an instance of a class.
        """
        # Checking solver
        self._GetPhysicsSolver().Check()
        self._GetAdjointSolver().Check()

        # Checking processes
        for process in self._GetListOfProcesses():
            process.Check()

    def RunSolutionLoop(self):
        self.InitializeTopologyOptimizationProblem()
        self.SolveTopologyOptimization()

    def InitializeTopologyOptimizationProblem(self):
        """
        This method Initializes the topology optimization problem solution
        """
        self.MpiPrint("\n------------------------------------------------------------------------------")
        self.MpiPrint(  "--| " + self.topology_optimization_name + " TOPOLOGY OPTIMIZATION PREPROCESSING")
        self.MpiPrint("------------------------------------------------------------------------------")
        self.MpiPrint("--|INITIALIZE|")
        self._GeometricalPreprocessing()
        self._InitializeOptimization()
        self._InitializeDomainDesign()
        self._PreprocessDerivatives() 

    def SolveTopologyOptimization(self):
        self.design_parameter_change = self.design_parameter_change_toll + 10
        end_solution = False
        self.first_iteration = True
        while (not end_solution):
            self.opt_it = self.opt_it+1
            self._GetSolver().main_model_part.ProcessInfo[KratosMultiphysics.STEP] = self.opt_it
            self.MpiPrint("\n------------------------------------------------------------------------------")
            self.MpiPrint("--| " + self.topology_optimization_name + " TOPOLOGY OPTIMIZATION SOLUTION LOOP. IT: " +  str(self.opt_it))
            self.MpiPrint("------------------------------------------------------------------------------")
            self.old_design_parameter = self.design_parameter
            self._SolveOptimizer()
            self._SolveTopologyOptimizationStepPhysics()
            self._EvaluateOptimizationProblem(print_results=True)
            end_solution = self._IsTopologyOptimizationSolutionEnd()
            if (end_solution):
                self.MpiPrint("\n------------------------------------------------------------------------------")
                self.MpiPrint("--| ENDING " + self.topology_optimization_name + " TOPOLOGY OPTIMIZATION SOLUTION LOOP")
                if (self.converged):
                    self.MpiPrint("--| ---> CONVERGED!")
                else:
                    self.MpiPrint("--| ---> Reached Max Number of Iterations")
                self._Remesh()
                self.MpiPrint("------------------------------------------------------------------------------\n")
            self._PrintOptimizationSolution()
            self.first_iteration = False

    def IsLastIteration(self):
        return (self.opt_it >= self.max_it)
            
    def _IsTopologyOptimizationSolutionEnd(self):
        design_parameter_converged = self._EvaluateDesignParameterChange()
        volume_constraint_valid = not (self.volume_constraint > 0.0)
        self.converged = design_parameter_converged and volume_constraint_valid
        # return not (((self.opt_it < self.max_it) and (not self.converged)) or (self.opt_it < self.min_it))
        return (self.opt_it >= self.max_it)

    def _SolveOptimizer(self):
        self.MpiPrint("\n--|OPTIMIZER|")
        self._SetTopologyOptimizationStage(3)
        if (self.first_iteration):
            self.MpiPrint("--|" + self.topology_optimization_stage_str + "| First Iteration: DO NOTHING")
            pass
        else:
            self.MpiPrint("--|" + self.topology_optimization_stage_str + "| SOLVE OPTIMIZER")
            opt_design_parameter = self._ExtractVariableInOptimizationDomain(self.design_parameter)
            self._SolveMMA(opt_design_parameter, self.n_opt_design_parameters, self.n_optimization_constraints, self.design_parameter_min_value, self.design_parameter_max_value, self.optimizer_max_outer_it, self.optimizer_kkt_tolerance)
    
    def _SolveTopologyOptimizationStepPhysics(self): 
        """
        This method organizs the physics solution of an optimization step
        """     
        ## Current Optimization Step Solutions 
        self._InitializeTopologyOptimizationStepPhysicsSolution()
        self._CheckMaterialProperties()
        self._UpdateRelevantPhysicsVariables()
        self._SolvePhysicsProblem() # NAVIER-STOKES PROBLEM SOLUTION
        self._PrintOnlyPhysics()
        if (self.first_iteration):
            self._SetFunctionalWeights()
            self._ResetFunctionalOutput()
        self._UpdateRelevantAdjointVariables()
        self._SolveAdjointProblem() # ADJOINT NAVIER-STOKES PROBLEM SOLUTION

    def _PrintOnlyPhysics(self, print=False):
        if (print):
            self.OptimizationOutputSolutionStep()
    
    def _GeometricalPreprocessing(self):
        """
        This method preprocess all the useful values and quantities for the topology optimization solution
        """
        self.MpiPrint("--|" + self.topology_optimization_stage_str + "| GEOMETRICAL PREPROCESSING")
        mp = self._GetComputingModelPart()
        self.dim = mp.ProcessInfo.GetValue(KratosMultiphysics.DOMAIN_SIZE)
        self.nodes_in_element = self.GetNumberOfNodesInElements()
        self._CreateNodesIdsDictionary()
        self._CreateElementsIdsDictionary()
        self._InitializeDomainSymmetry()
        self._ComputeDomainSize()

    def GetNumberOfNodesInElements(self):
        num_nodes_elements = 0
        if (len(self._GetMainModelPart().Elements) > 0):
            for elem in self._GetMainModelPart().Elements:
                num_nodes_elements = len(elem.GetNodes())
                break
        return num_nodes_elements
    
    def _InitializeOptimization(self):  
        self.MpiPrint("--|" + self.topology_optimization_stage_str + "| INITIALIZE OPTIMIZATION")
        self.SetTimeSolutionStorage()
        self._ComputeDesignParameterFilterUtilities()
        self._InitializeOptimizationSolutionStabilization()
        self.opt_it = 0
        self._SetDesignParameterChangeTolerance()
        self.n_optimization_constraints = 0  
        self._InitializeOptimizerSettings()
        self._InitializeConstraints()
        self._InitializeFunctionalEvaluation()
        self._InitializeRemeshing()

    def SetTimeSolutionStorage(self):
        # velocity
        self.velocity_solutions         = np.zeros((self.n_time_steps, self.n_nodes, self.dim))
        self.adjoint_velocity_solutions = np.zeros((self.n_time_steps, self.n_nodes, self.dim))
        # pressure
        self.pressure_solutions         = np.zeros((self.n_time_steps, self.n_nodes))
        self.adjoint_pressure_solutions = np.zeros((self.n_time_steps, self.n_nodes))
    
    def _SetDesignParameterChangeTolerance(self):
        self.design_parameter_change_toll = self.optimization_settings["change_tolerance"].GetDouble()

    def _InitializeOptimizerSettings(self):
        self.optimizer_settings = self.optimization_settings["optimizer_settings"]
        self.optimizer_name = self.optimizer_settings["optimizer"].GetString()
        self.design_parameter_min_value = self.optimizer_settings["values_range"][0].GetDouble()
        self.design_parameter_max_value = self.optimizer_settings["values_range"][1].GetDouble()
        self.optimizer_max_outer_it = self.optimizer_settings["max_outer_it"].GetInt()
        self.optimizer_kkt_tolerance = self.optimizer_settings["kkt_tolerance"].GetDouble()

    def _InitializeConstraints(self):
        self.constraints_settings = self.optimization_settings["optimization_problem_settings"]["constraints_settings"]
        self._InitializeVolumeConstraint()
        self._InitializeOtherConstraints()
        self._ResetConstraints()

    def _InitializeVolumeConstraint(self):
        self.volume_constraint_settings = self.constraints_settings["volume_constraint_settings"]
        self.volume_constraint_id = self.n_optimization_constraints
        self.n_optimization_constraints += 1
        if (self.volume_constraint_settings["fluid_or_solid"].GetString().lower() == "solid"):
            self.is_fluid_volume_constraint = False
        else:
            self.is_fluid_volume_constraint = True
        self._SetMaxDomainVolumeFraction()

    def _InitializeOtherConstraints(self):
        # Set use the specific contraint to 'False' for all the custom implemented constraints
        self.use_wss_constraint = False

        self.use_other_constraints = self.constraints_settings["use_other_constraints"].GetBool()
        if (self.use_other_constraints):
            self._InitializeWSSConstraint()

    def _InitializeWSSConstraint(self):
        self.wss_constraint_settings = self.constraints_settings["other_constraints_list"]["WSS_constraint_settings"]
        self.use_wss_constraint = self.wss_constraint_settings["use_WSS_constraint"].GetBool()
        if (self.use_wss_constraint):
            self.wss_constraint_id = self.n_optimization_constraints
            self.n_optimization_constraints += 1
            self.min_wss = self.wss_constraint_settings["min_WSS"].GetDouble()

    def _ResetConstraints(self):
        self.constraints = np.zeros((self.n_optimization_constraints,1))  
        self.constraints_derivatives_wrt_design = np.zeros((self.n_optimization_constraints, self.n_opt_design_parameters))

    def _InitializeFunctionalEvaluation(self):
        self.functional_weights_imported = False
        self._InitializeFunctionalsInTime()
        self._ImportFunctionalWeights()

    def _InitializeFunctionalsInTime(self):
        self._InitializeFluidFunctionalsInTime()
    
    def _InitializeFluidFunctionalsInTime(self):
        self._InitializeResistanceFunctionalInTime()
        self._InitializeStrainRateFunctionalInTime()
        self._InitializeVorticityFunctionalsInTime()

    def _InitializeResistanceFunctionalInTime(self):
        self.resistance_functionals_in_delta_time  = np.zeros(self.n_time_steps)
        self.resistance_functional_time_steps_integration_weights, self.resistance_functional_time_info = self._InitializeFunctionalTimeIntegrationWeights(functional_physics="fluid_functionals", functional_name="resistance")

    def _InitializeStrainRateFunctionalInTime(self):
        self.strain_rate_functionals_in_delta_time = np.zeros(self.n_time_steps)
        self.strain_rate_functional_time_steps_integration_weights, self.strain_rate_functional_time_info = self._InitializeFunctionalTimeIntegrationWeights(functional_physics="fluid_functionals", functional_name="strain_rate")
    
    def _InitializeVorticityFunctionalsInTime(self):
        self.vorticity_functionals_in_delta_time   = np.zeros(self.n_time_steps)
        self.vorticity_functional_time_steps_integration_weights, self.vorticity_functional_time_info = self._InitializeFunctionalTimeIntegrationWeights(functional_physics="fluid_functionals", functional_name="vorticity")

    def _InitializeFunctionalTimeIntegrationWeights(self, functional_physics, functional_name):
        functional_time_steps_integration_weights = np.zeros(self.n_time_steps)
        times = {}
        times["times"] = np.asarray([self.start_time, self.end_time])
        times["start_span_times"] = np.asarray([self.start_time, self.end_time])
        times["start_span_time_steps_ids"] = [0,1]
        times["end_span_times"] = np.asarray([self.start_time, self.end_time])
        times["end_span_time_steps_ids"] = [0,1]
        if self.IsUnsteadySolution():
            [start_time, end_time] = self._GetFunctionalEvaluationTimeInterval(functional_physics, functional_name)
            times["times"] = np.asarray([start_time, end_time])
            dt = end_time-start_time    
            start_id = np.searchsorted(self.simulation_times, start_time, side="right")
            end_id = np.searchsorted(self.simulation_times, end_time, side="right")
            if (end_id == self.n_time_steps+1):
                end_id = end_id-1
            if end_id > start_id:
                start_span_times_ids = [start_id-1, start_id]  
                end_span_times_ids   = [end_id-1, end_id]  
                start_span_time_steps_ids = [id-1 for id in start_span_times_ids] # time_step ids that multiply the the start_span_times in time integral evaluations (functional at time_id=n multiplies the time_step_id=n-1)
                times["start_span_time_steps_ids"] = start_span_time_steps_ids
                end_span_time_steps_ids   = [id-1 for id in end_span_times_ids] # time_step ids that multiply the the end_span_times in time integral evaluations
                times["end_span_time_steps_ids"] = end_span_time_steps_ids
                start_span_times = [self.simulation_times[start_span_times_ids[0]], self.simulation_times[start_span_times_ids[1]]] # interval containing time_interval[0]
                times["start_span_times"] = np.asarray(start_span_times)
                end_span_times   = [self.simulation_times[end_span_times_ids[0]], self.simulation_times[end_span_times_ids[1]]] # interval containing time_interval[1]
                times["end_span_times"] = np.asarray(end_span_times)
                dt_start_span = start_span_times[1]-start_span_times[0]
                dt_end_span   = end_span_times[1]-end_span_times[0] 
                start_eval_weights = [(start_span_times[1]-start_time)/dt_start_span, (start_time-start_span_times[0])/dt_start_span]
                end_eval_weights   = [(end_span_times[1]-end_time)/dt_end_span, (end_time-end_span_times[0])/dt_end_span]
                if (start_span_time_steps_ids[0] >= 0):
                    # check to store the first start_span_time_step_weight only if it is not relative to the starting time, 
                    # at the starting time the solution is assumed to be zero and thererefore is not integrated
                    functional_time_steps_integration_weights[start_span_time_steps_ids[0]] += (start_eval_weights[0]**2)  * dt_start_span / 2.0
                functional_time_steps_integration_weights[start_span_time_steps_ids[1]]     += (1.0+start_eval_weights[1]) * start_eval_weights[0] * dt_start_span / 2.0
                current_id = start_span_time_steps_ids[1]+1
                while current_id <= end_span_time_steps_ids[0]:
                    functional_time_steps_integration_weights[current_id-1]   += self.time_steps[current_id] / 2.0
                    functional_time_steps_integration_weights[current_id]     += self.time_steps[current_id] / 2.0
                    current_id += 1
                functional_time_steps_integration_weights[end_span_time_steps_ids[0]]   += (1.0+end_eval_weights[0]) * end_eval_weights[1] * dt_end_span / 2.0
                functional_time_steps_integration_weights[end_span_time_steps_ids[1]]   += (end_eval_weights[1]**2)  * dt_end_span / 2.0
            else: # end_id = start_id
                span_times_ids = [start_id-1, start_id]  # unique span ids array
                span_time_steps_ids = [id-1 for id in span_times_ids]
                times["start_span_time_steps_ids"] = span_time_steps_ids
                times["end_span_time_steps_ids"] = span_time_steps_ids
                span_times = [self.simulation_times[span_times_ids[0]], self.simulation_times[span_times_ids[1]]]
                times["start_span_times"] = np.asarray(span_times)
                times["end_span_times"] = np.asarray(span_times)
                dt_span = span_times[1]-span_times[0]
                start_eval_weights = [(span_times[1]-start_time)/dt_span, (start_time-span_times[0])/dt_span]
                end_eval_weights   = [(span_times[1]-end_time)/dt_span, (end_time-span_times[0])/dt_span]
                if (span_time_steps_ids[0] >= 0):
                    # check to store the first start_span_time_step_weight only if it is not relative to the starting time, 
                    # at the starting time the solution is assumed to be zero and thererefore is not integrated
                    functional_time_steps_integration_weights[span_time_steps_ids[0]] += (start_eval_weights[0]+end_eval_weights[0]) * dt / 2.0
                functional_time_steps_integration_weights[span_time_steps_ids[1]] += (start_eval_weights[1]+end_eval_weights[1]) * dt / 2.0
        else: # steady simulation, no time integration
            functional_time_steps_integration_weights = [1.0]
        return functional_time_steps_integration_weights, times

    def _GetFunctionalEvaluationTimeInterval(self, functional_physics, functional_name):
        physics_time_interval_list = self.optimization_settings["optimization_problem_settings"]["functional_weights"][functional_physics]["time_interval"]
        functional_time_interval_list = self.optimization_settings["optimization_problem_settings"]["functional_weights"][functional_physics][functional_name]["time_interval"]
        physics_start_time, physics_end_time = self._GetTimeIntervalFromSettingsList(physics_time_interval_list, functional_err_msg=functional_physics)
        functional_start_time, functional_end_time = self._GetTimeIntervalFromSettingsList(functional_time_interval_list, functional_err_msg=functional_physics+"->"+functional_name)
        if physics_start_time > functional_end_time:
            info_msg = "Calling '_GetFunctionalEvaluationTimeInterval' for functional: " + functional_physics + "->" + functional_name + " setting physics_start_time > functional_end_time"
            raise RuntimeError("FluidTransportTopologyOptimizationAnalysis: " + info_msg)
        elif functional_start_time > physics_end_time:
            info_msg = "Calling '_GetFunctionalEvaluationTimeInterval' for functional: " + functional_physics + "->" + functional_name + " setting functional_start_time > physics_end_time"
        start_time_value = max(physics_start_time, functional_start_time)
        end_time_value = min(physics_end_time, functional_end_time)
        if start_time_value > end_time_value:
            info_msg = "Calling '_GetFunctionalEvaluationTimeInterval' for functional: " + functional_physics + "->" + functional_name + " setting start_time_value > end_time_value"
            raise RuntimeError("FluidTransportTopologyOptimizationAnalysis: " + info_msg)
        return start_time_value, end_time_value
    
    def _GetTimeIntervalFromSettingsList(self, interval_list, functional_err_msg="Not Specified"):
        if interval_list[0].IsString():
            start_time_value = interval_list[0].GetString()
            if start_time_value == "Start":
                start_time_value = self.start_time
            else:
                info_msg = "Calling '_GetTimeIntervalFromSettings' for functional: " + functional_err_msg + " setting a start time using a string_type != 'String'"
                raise RuntimeError("FluidTransportTopologyOptimizationAnalysis: " + info_msg)
        else:
            start_time_value = max(interval_list[0].GetDouble(), self.start_time)
        if interval_list[1].IsString():
            end_time_value = interval_list[1].GetString()
            if end_time_value == "End":
                end_time_value = self.end_time
            else:
                info_msg = "Calling '_GetTimeIntervalFromSettings' for functional: " + functional_err_msg + " setting a end time using a string_type != 'End'"
                raise RuntimeError("FluidTransportTopologyOptimizationAnalysis: " + info_msg)
        else:
            end_time_value = min(interval_list[1].GetDouble(), self.end_time)
        if start_time_value > end_time_value:
            info_msg = "Calling '_GetTimeIntervalFromSettings' for functional: " + functional_err_msg + " setting start_time > end_time"
            raise RuntimeError("FluidTransportTopologyOptimizationAnalysis: " + info_msg)
        return start_time_value, end_time_value

    def _InitializeDomainDesign(self):
        self.MpiPrint("--|" + self.topology_optimization_stage_str + "| INITIALIZE DOMAIN DESIGN")
        self._InitializeDomainDesignParameter()
        self._InitializePhysicsParameters()

    def _InitializePhysicsParameters(self):
        self._InitializeResistance()
    
    def _SetMaxDomainVolumeFraction(self):
        self.max_volume_fraction = self.volume_constraint_settings["max_volume_fraction"].GetDouble()

    def _InitializeDomainDesignParameter(self):
        """
        This method will handle the design parameter initialization across the whole domain.
        self.design_parameter is defined in every node of the mesh.
        The optimization process will update only the value of the nodes belonging to the optimization_domain sub model part
        """
        mask = self._GetOptimizationDomainNodesMask()
        # initial_design_settings = self.optimization_settings["optimization_problem_settings"]["initial_design_settings"]
        # default_initial_design_value = initial_design_settings["default_initial_design"].GetDouble()
        self.volume_fraction = self.optimization_domain_initial_value
        self.design_parameter = np.zeros(self.n_nodes) + self.non_optimization_domain_initial_value
        self.design_parameter[mask] = self.optimization_domain_initial_value
        self._SetDesignParameterCustomInitialDesign()
        counter = 0
        for node in self._GetLocalMeshNodes():
            design = self.design_parameter[counter]
            counter += 1
            node.SetSolutionStepValue(KratosMultiphysics.DESIGN_PARAMETER, design)
            distance = design-self.remeshing_levelset
            node.SetSolutionStepValue(KratosMultiphysics.DISTANCE, distance)

    def _SetDesignParameterCustomInitialDesign(self):
        if self.use_custom_initial_design:
            custom_initial_designs_list = self.custom_initial_design_settings["custom_initial_design"]
            for custom_initial_design in custom_initial_designs_list:
                model_part_nodes_ids = self._ExtractListOfNodesFromNodesDictionary(self._GetSubModelPart(self._GetMainModelPart(), custom_initial_design["model_part_name"].GetString()))
                initial_design_value = custom_initial_design["initial_value"].GetDouble()
                self.design_parameter[model_part_nodes_ids] = initial_design_value
        else:
            pass

    def _InitializeDomainSymmetry(self):
        self.symmetry_settings = self.optimization_settings["symmetry_settings"]
        self.symmetry_enabled = self.symmetry_settings["symmetry"].GetBool()
        self.symmetry_model_part_name = self.symmetry_settings["model_part_name"].GetString()
        if (self.symmetry_enabled):
            sym_mp = self._GetSubModelPart(self._GetMainModelPart(), self.symmetry_model_part_name)
            if (sym_mp is None):
                raise RuntimeError("Trying to exploit symmetry in a non existent symmetry plane model part")
            else:
                self.symmetry_nodes_mask = np.zeros(len(self._GetLocalMeshNodes(sym_mp)), dtype=int)
                count = 0
                for node in self._GetLocalMeshNodes(sym_mp):
                    self.symmetry_nodes_mask[count] = self.nodes_ids_global_to_local_partition_dictionary[node.Id]
                    count +=1
        else:
            self.symmetry_nodes_mask = np.asarray([])

    def _CorrectNodalDomainSizeWithSymmetry(self):
        if (self.symmetry_enabled):
            self.nodal_domain_sizes[self.symmetry_nodes_mask] *= 2.0

    def _ComputeScalarVariableNodalGradient(self, scalar_variable, gradient_variable):
        mp = self._GetComputingModelPart()
        gradient_process = ComputeNodalGradientProcess(mp, scalar_variable, gradient_variable, KratosMultiphysics.NODAL_AREA)
        gradient_process.Execute()

    def _ComputeDesignParameterFilterUtilities(self):
        self._ComputeDesignParameterDiffusiveFilterUtilities()
        self._ComputeDesignParameterProjectiveFilterUtilities()

    def _InitializeOptimizationSolutionStabilization(self):
        self.optimization_solution_stabilization_settings = self.optimization_settings["solution_stabilization_settings"]
        self._InitializeAdjointPhysicsSolutionStabilization()

    def _InitializeAdjointPhysicsSolutionStabilization(self):
        self._InitializeAdjointFluidSolutionStabilization()

    def _InitializeAdjointFluidSolutionStabilization(self):
        self._InitializeAdjointSolutionViscosityAdaptation()

    def _InitializeAdjointSolutionViscosityAdaptation(self):
        self.adjoint_viscosity_adaptation_settings = self.optimization_solution_stabilization_settings["adjoint_viscosity_adaptation_settings"]
        self.use_adjoint_viscosity_adaptation = self.adjoint_viscosity_adaptation_settings["use_adjoint_viscosity_adaptation"].GetBool()
        if self.use_adjoint_viscosity_adaptation:
            self.adjoint_viscosity_adaptation_factor = self.adjoint_viscosity_adaptation_settings["adjoint_viscosity_adaptation_factor"].GetDouble()
            if (abs(self.adjoint_viscosity_adaptation_factor-1.0) > 1e-10):
                self.physics_viscosity = self._GetViscosity()
                self.adjoint_viscosity = self.physics_viscosity * self.adjoint_viscosity_adaptation_factor
            else:
                self.use_adjoint_viscosity_adaptation = False

    def _InitializeOptimizationDomainSettings(self):
        """
        Reads and initializes the settings related to the optimization domain.
        It retrieves:
        - The model part names and initial values for both the optimization and non-optimization domains,
        ensuring that the initial values remain within [0.0, 1.0].
        - The configuration for a possible custom initial design, activated through `use_custom_initial_design`.
        """
        optimization_domain_settings = self.optimization_settings["optimization_domain_settings"]
        self.optimization_domain_name = optimization_domain_settings["optimization_domain"]["model_part_name"].GetString()
        self.optimization_domain_initial_value = max(0.0, min(1.0, optimization_domain_settings["optimization_domain"]["initial_value"].GetDouble()))
        self.non_optimization_domain_name = optimization_domain_settings["non_optimization_domain"]["model_part_name"].GetString()
        self.non_optimization_domain_initial_value = max(0.0, min(1.0, optimization_domain_settings["non_optimization_domain"]["initial_value"].GetDouble()))
        self.custom_initial_design_settings = optimization_domain_settings["custom_initial_design_settings"]
        self.use_custom_initial_design = self.custom_initial_design_settings["use_custom_initial_design"].GetBool()

    def _CorrectNodalOptimizationDomainSizeWithSymmetry(self):
        if (self.symmetry_enabled):
            self.nodal_optimization_domain_sizes[self.symmetry_nodes_mask] *= 2.0

    def _GetModelPartNodesSubset(self, model_part, node_ids):
        """
        Returns a list of nodes from the given model part corresponding to the specified node IDs.
        Each node is retrieved via `model_part.GetNode(node_id)`.
        The input `node_ids` must be global node IDs (Kratos uses global IDs also in MPI runs).
        In non-MPI runs, global and local IDs coincide.
        """
        return [model_part.GetNode(node_id) for node_id in node_ids]
    
    def _GetModelPartNodesIds(self, model_part):
        """
        Returns the list of global node IDs belonging to the given model part.
        The nodes are obtained from `_GetLocalMeshNodes(model_part)`, which handles both
        MPI and non-MPI cases by returning the locally owned nodes in parallel runs
        or all nodes in serial runs.
        """
        model_part_nodes = self._GetLocalMeshNodes(model_part)
        return [node.Id for node in model_part_nodes]
            
    def _ExtractListOfNodesFromNodesDictionary(self, model_part):
        """
        Returns the list of local node indices corresponding to the nodes in the given model part.
        Global node IDs from the model part are converted to local indices using the
        `nodes_ids_global_to_local_partition_dictionary` created in `_CreateNodesIdsDictionary`.
        If the simulation is not running in MPI, the returned list simply matches the global node IDs,
        since the mapping is one-to-one.
        """
        model_part_nodes_ids = self._GetModelPartNodesIds(model_part)
        nodes_list = [self.nodes_ids_global_to_local_partition_dictionary[node_id] for node_id in model_part_nodes_ids]
        return nodes_list

    def _ExtractVariableInOptimizationDomain(self, variable):
        return variable[self.optimization_domain_nodes_mask]
    
    def _InsertDesignParameterFromOptimizationDomain(self, optimization_domain_design_parameter):
        mask = self._GetOptimizationDomainNodesMask()
        design_parameter = self.design_parameter.copy()
        design_parameter[mask] = optimization_domain_design_parameter
        return design_parameter
    
    def _PreprocessDerivatives(self):
        """
        This method preprocess the quantities for easy derivatives evaluation
        """
        self.MpiPrint("--|" + self.topology_optimization_stage_str + "| PREPROCESS DERIVATIVES")
        self._PreprocessGradient()

    def _PreprocessGradient(self):
        """Avoided method in mpi: the velocity gradient is evaluated on the nodes insted of elements, maybe this can be improved"""
        pass

    def _GetShapeFunctionsDerivatives(self, element):
        gradients = np.asarray(element.GetGeometry().ShapeFunctionDerivatives(1,0))
        jacobian  = np.asarray(element.GetGeometry().Jacobian(0))
        inv_jacobian = np.linalg.inv(jacobian)
        global_gradient = gradients @ inv_jacobian
        return global_gradient 
    
    def _InitializeTopologyOptimizationStepPhysicsSolution(self):
        """
        This method Re-Initializes the physical values in order to be consistent with the solution of the current optimization step
        """
        self.MpiPrint("\n--|INITIALIZE OPTIMIZATION STEP PHYSICS SOLUTION|")
        self._SetTopologyOptimizationStage(0)
        if (not self.first_iteration):
                self._ReInitializePhysics()
        self._UpdateDesignParameterAndPhysicsParameters(self.design_parameter)

    def _UpdateDesignParameterAndPhysicsParameters(self, design_parameter):
        self._UpdateDesignParameter(design_parameter)
        self._UpdatePhysicsParameters()
    
    def _UpdatePhysicsParameters(self):
        self._UpdateResistance()

    def UpdatePhysicsParametersVariables(self):
        self._UpdateResistanceVariable()

    def _UpdateDesignParameter(self, design_parameter):
        self.MpiPrint("--|" + self.topology_optimization_stage_str + "| UPDATE DESIGN PARAMETER")
        self.design_parameter_base = design_parameter
        self._ApplyDesignParameterDiffusiveFilter(design_parameter)
        self._ApplyDesignParameterProjectiveFilter(self.design_parameter_filtered)
        self.design_parameter = self.design_parameter_projected
        self._UpdateDesignParameterVariable()

    def _UpdateDesignParameterVariable(self):
        for node in self._GetLocalMeshNodes():
            design = self.design_parameter[self.nodes_ids_global_to_local_partition_dictionary[node.Id]]
            node.SetSolutionStepValue(KratosMultiphysics.DESIGN_PARAMETER, design)
            distance = design-self.remeshing_levelset
            node.SetSolutionStepValue(KratosMultiphysics.DISTANCE, distance)

    def _ResetPhysicsParameters(self):
        self._ResetResistance()

    def _InitializeResistance(self):
        self.MpiPrint("--|" + self.topology_optimization_stage_str + "| INITIALIZE RESISTANCE")
        self.resistance_parameters = self.physics_parameters_settings["resistance"] 
        self._ResetResistance()
        self._UpdateResistanceVariable()

    def _ResetResistance(self):
        self.resistance = np.zeros(self.n_nodes)
        self.resistance_derivative_wrt_design = np.zeros(self.n_nodes)
    
    def _UpdateResistance(self):
        """
        This method handles the resistance update.
        """
        self.MpiPrint("--|" + self.topology_optimization_stage_str + "| UPDATE RESISTANCE")
        self.resistance, self.resistance_derivative_wrt_design_base = self._ComputeResistance(self.design_parameter)
        self._UpdateResistanceDesignDerivative()
        self._UpdateResistanceVariable()

    def _ComputeResistance(self, design_parameter):
        self._EvaluatePhysicsParameterCurrentValue(self.resistance_parameters)
        return self._ComputePhysicsParameter(self.resistance_parameters, design_parameter)
    
    def _ComputePhysicsParameter(self, physics_parameters, design_parameter):
        parameter = np.zeros(self.n_nodes)
        parameter_derivative_wrt_design_base = np.zeros(self.n_nodes)
        for subdomain_physics_parameters in physics_parameters.values():
            domain = subdomain_physics_parameters["domain"].GetString()
            mp = self._GetSubModelPart(self._GetMainModelPart(), domain)
            if (mp is not None):
                nodes_ids = self._ExtractListOfNodesFromNodesDictionary(mp)
            else:
                nodes_ids = np.arange(self.n_nodes)
            interpolation_method = (subdomain_physics_parameters["interpolation_method"].GetString()).lower()
            if (interpolation_method == "hyperbolic"):
                computed_parameter, computed_parameter_derivative_wrt_design_base = self._ComputeHyperbolicPhysicsParameter(subdomain_physics_parameters, design_parameter)
            elif (interpolation_method == "polynomial"):
                computed_parameter, computed_parameter_derivative_wrt_design_base = self._ComputePolynomialPhysicsParameter(subdomain_physics_parameters, design_parameter)
            else:
                raise RuntimeError("WRONG PHYSICS PARAMETER INTERPOLATION METHOD", "Running '_ComputePhysicsParameter' with the wrong interpolation method.")
            parameter[nodes_ids] = computed_parameter[nodes_ids]
            parameter_derivative_wrt_design_base[nodes_ids] = computed_parameter_derivative_wrt_design_base[nodes_ids]
        return parameter, parameter_derivative_wrt_design_base
        
    def _ComputeHyperbolicPhysicsParameter(self, physics_parameters, design_parameter):
        value_void = physics_parameters["value_void"].GetDouble()
        value_full = physics_parameters["value_full"].GetDouble()
        slope_initial, slope_final, start_it, end_it = self._GetPhysicsParameterSlopeScaling(physics_parameters["interpolation_slope"])
        if (slope_initial > 1.0):
            slope_initial = 1.0
        if (slope_final > slope_initial):
            slope_final = slope_initial
        eff_slope = self._LinearInterpolationOverIterations(slope_initial, slope_final, start_it, end_it)
        if (value_void <= value_full): 
            physics_parameter = value_void + (value_full-value_void)*(eff_slope*design_parameter)/(eff_slope+1-design_parameter)
            physics_parameter_derivative_wrt_design_base = (value_full-value_void)*(eff_slope*(eff_slope+1))/((eff_slope+1-design_parameter)**2)
        else:
            physics_parameter = value_full - (value_full-value_void)*(eff_slope*(1-design_parameter))/(eff_slope+design_parameter)
            physics_parameter_derivative_wrt_design_base = (value_full-value_void)*(eff_slope*(eff_slope+1))/((eff_slope+design_parameter)**2)
        return physics_parameter, physics_parameter_derivative_wrt_design_base
    
    def _ComputePolynomialPhysicsParameter(self, physics_parameters, design_parameter):
        value_void = physics_parameters["value_void"].GetDouble()
        value_full = physics_parameters["value_full"].GetDouble()
        power_initial, power_final, start_it, end_it = self._GetPhysicsParameterSlopeScaling(physics_parameters["interpolation_slope"])
        if (power_initial < 1.0):
            power_initial = 1.0
        if (power_final < power_initial):
            power_final = power_initial
        eff_power = self._LinearInterpolationOverIterations(power_initial, power_final, start_it, end_it)
        if (value_void <= value_full): 
            physics_parameter = value_void + (value_full-value_void)*(design_parameter**eff_power)
            physics_parameter_derivative_wrt_design_base = eff_power*(value_full-value_void)*(design_parameter**(eff_power-1))
        else:
            physics_parameter = value_full - (value_full-value_void)*((1.0-design_parameter)**eff_power)
            physics_parameter_derivative_wrt_design_base = eff_power*(value_full-value_void)*((1.0-design_parameter)**(eff_power-1))
        return physics_parameter, physics_parameter_derivative_wrt_design_base
    
    def _LinearInterpolationOverIterations(self, initial_value, final_value, start_it, end_it):
        if (self.opt_it < start_it):
            return initial_value
        elif (self.opt_it < end_it):
            return initial_value + (final_value-initial_value)*(self.opt_it-start_it)/(end_it-start_it)
        else:
            return final_value
        
    def _GetPhysicsParameterSlopeScaling(self, interpolation_slope_parameters):
        slope_initial = interpolation_slope_parameters["initial_slope"].GetDouble()
        slope_final   = interpolation_slope_parameters["final_slope"].GetDouble()
        start_it = interpolation_slope_parameters["iterations"][0].GetInt()
        end_it   = interpolation_slope_parameters["iterations"][1].GetInt()
        return slope_initial, slope_final, start_it, end_it
    
    def _EvaluatePhysicsParameterCurrentValue(self, physics_parameters):
        for subdomain_physics_parameters in physics_parameters.values():
            # Adjust void value
            if (subdomain_physics_parameters["change_value_void"]["change_value"].GetBool()):
                start_it      = subdomain_physics_parameters["change_value_void"]["iterations"][0].GetInt()
                end_it      = subdomain_physics_parameters["change_value_void"]["iterations"][1].GetInt()
                initial_value = subdomain_physics_parameters["change_value_void"]["initial_value"].GetDouble()
                final_value = subdomain_physics_parameters["change_value_void"]["final_value"].GetDouble()
                if (self.opt_it < self.min_it):
                    current_value = initial_value
                elif (self.opt_it < end_it):
                    current_value = initial_value + (final_value-initial_value)*(self.opt_it-start_it)/(end_it-start_it)
                else:
                    current_value = final_value
                subdomain_physics_parameters["value_void"].SetDouble(current_value)
            # Adjust full value
            if (subdomain_physics_parameters["change_value_full"]["change_value"].GetBool()):
                start_it      = subdomain_physics_parameters["change_value_full"]["iterations"][0].GetInt()
                end_it      = subdomain_physics_parameters["change_value_full"]["iterations"][1].GetInt()
                initial_value = subdomain_physics_parameters["change_value_full"]["initial_value"].GetDouble()
                final_value = subdomain_physics_parameters["change_value_full"]["final_value"].GetDouble()
                if (self.opt_it < start_it):
                    current_value = initial_value
                elif (self.opt_it < end_it):
                    current_value = initial_value + (final_value-initial_value)*(self.opt_it-start_it)/(end_it-start_it)
                else:
                    current_value = final_value
                subdomain_physics_parameters["value_full"].SetDouble(current_value)

    def _UpdateResistanceDesignDerivative(self):
        resistance_derivative_wrt_design_projected = self.resistance_derivative_wrt_design_base * self.design_parameter_projected_derivatives
        self.resistance_derivative_wrt_design = resistance_derivative_wrt_design_projected
    
    def _UpdateResistanceVariable(self):
        self._GetSolver()._UpdateResistanceVariable(self.resistance)

    def _SolvePhysicsProblem(self):  
        """
        This method executes a single NS solution of the Topology Optimization problem loop
        """
        self._SetTopologyOptimizationStage(1)# PHYSICS PROBLEM SOLUTION 
        self._RunStageSolutionLoop()  
    
    def _SolveAdjointProblem(self):  
        """
        This method executes a single ADJ_NS solution of the Topology Optimization problem loop
        """
        self._SetTopologyOptimizationStage(2) # ADJOINT PROBLEM SOLUTION
        self._InitializeAdjointPhysicsSolutionLoopStabilization()
        self._RunStageSolutionLoop() 
        self._FinalizeAdjointPhysicsSolutionLoopStabilization()

    def _RunStageSolutionLoop(self):
        """
        This method executes a single physics solution of the Topology Optimization problem loop
        N.B.: must be called after the creation of the fluid-transport model part
        problem_stage = 1: Navier Stokes solution
        problem_stage = 2: Adjoint Navier-Stokes solution
        """
        problem_stage = self.GetTopologyOptimizationStage()
        if ((problem_stage == 1)): 
            self.MpiPrint("\n--|PHYSICS SOLUTION|")
        elif (problem_stage == 2):
            self.MpiPrint("\n--|ADJOINT SOLUTION|")
        else:   
            self.MpiPrint("--| UNKNOWN SOLUTION |")
        top_opt_stage_str = self.topology_optimization_stage_str
        self.SetupStageSolutionLoopTimeSettings()
        self.MpiPrint("--|" + top_opt_stage_str + "| START SOLUTION LOOP")
        while self.KeepAdvancingSolutionLoop():
            self.MpiPrint("--|" + top_opt_stage_str + "| ADVANCE TIME STEP")
            self.time = self._AdvanceTime()
            self.MpiPrint("--|" + top_opt_stage_str + "| INITIALIZE SOLUTION TIME STEP")
            self.InitializeSolutionStep()
            self.MpiPrint("--|" + top_opt_stage_str + "| UPDATE PHYSICS PARAMETERS TIME STEP")
            self.UpdatePhysicsParametersVariablesAndSynchronize()
            self.MpiPrint("--|" + top_opt_stage_str + "| PREDICT")
            self._GetSolver().Predict()
            self.MpiPrint("--|" + top_opt_stage_str + "| SOLVE SOLUTION TIME STEP")
            is_converged = self._GetSolver().SolveSolutionStep()
            self.MpiPrint("--|" + top_opt_stage_str + "| CHECK CONVERGENCE")
            self.__CheckIfSolveSolutionStepReturnsAValue(is_converged)
            self.MpiPrint("--|" + top_opt_stage_str + "| FINALIZE SOLUTION TIME STEP")
            self.FinalizeSolutionStep()
        self.MpiPrint("--|" + top_opt_stage_str + "| END SOLUTION LOOP")

    def SetupStageSolutionLoopTimeSettings(self):
        self.time = self.start_time
        self._GetSolver()._SetStartingTime(self.start_time)
        self.time_step_counter = 0
        self.ResetPhysicsTimeStepVariables()
        
    def ResetPhysicsTimeStepVariables(self):
        self._ResetFluidTimeStepVariables()
    
    def _ResetFluidTimeStepVariables(self):
        if (self.IsPhysicsStage()):
            self._GetPhysicsSolver().main_model_part.ProcessInfo[KratosCFD.FLUID_TOP_OPT_NS_STEP] = 0
        elif (self.IsAdjointStage()):
            self._GetAdjointSolver().main_model_part.ProcessInfo[KratosCFD.FLUID_TOP_OPT_ADJ_NS_STEP] = 0

    def _InitializeAdjointPhysicsSolutionLoopStabilization(self):
        self._InitializeAdjointFluidSolutionLoopStabilization()

    def _InitializeAdjointFluidSolutionLoopStabilization(self):
        self._InitializeAdjointSolutionLoopViscosityAdaptation()
    
    def _InitializeAdjointSolutionLoopViscosityAdaptation(self):
        if (self.use_adjoint_viscosity_adaptation):
            for elem in self._GetMainModelPart().GetCommunicator().LocalMesh().Elements:
                elem.Properties[KratosMultiphysics.DYNAMIC_VISCOSITY] = self.adjoint_viscosity

    def _FinalizeAdjointPhysicsSolutionLoopStabilization(self):
        self._FinalizeAdjointFluidSolutionLoopStabilization()

    def _FinalizeAdjointFluidSolutionLoopStabilization(self):
        self._FinalizeAdjointSolutionLoopViscosityAdaptation()

    def _FinalizeAdjointSolutionLoopViscosityAdaptation(self):
        if (self.use_adjoint_viscosity_adaptation):
            for elem in self._GetMainModelPart().GetCommunicator().LocalMesh().Elements:
                elem.Properties[KratosMultiphysics.DYNAMIC_VISCOSITY] = self.physics_viscosity

    def KeepAdvancingSolutionLoop(self):
        """This method specifies the stopping criteria for breaking the solution loop
        It can be overridden by derived classes
        """
        state =  (abs(self.end_time-self.time) > 1e-9*self.end_time)
        if (self.time_step_counter > self.n_time_steps):
            raise RuntimeError("Advancing in Solution Loop with: (self.time_step_counter > self.n_time_steps")
        return state

    def _AdvanceTime(self):
        """
        This methods hadls the time for the Topology Optimization problems, currently its purpos it's just to make thinks work
        """
        self.time_step_counter += 1
        if self.IsUnsteadySolution():
            time = super()._AdvanceTime()
        else:
            time = self.start_time + self.time_step_base
        return time
    
    def InitializeSolutionStep(self):
        super().InitializeSolutionStep()
        if self.IsPhysicsStage():
            self.InitializePhysicsSolutionStep()
        if self.IsAdjointStage():
            self.UpdateFunctionalWeights()
            self.InitializeAdjointPhysicsSolutionStep()
        for process in self._GetListOfTimeOutputProcesses():
            process.ExecuteInitializeSolutionStep()

    def UpdateFunctionalWeights(self):
        if self.IsUnsteadySolution():
            time_step_functional_weights = self.functional_weights_in_time[self.n_time_steps-self.time_step_counter,:] * self.functional_weights
            self._GetComputingModelPart().ProcessInfo.SetValue(KratosMultiphysics.FUNCTIONAL_WEIGHTS, time_step_functional_weights)
        else:
            self._GetComputingModelPart().ProcessInfo.SetValue(KratosMultiphysics.FUNCTIONAL_WEIGHTS, self.functional_weights)
        self._GetMainModelPart().GetCommunicator().SynchronizeNonHistoricalVariable(KratosMultiphysics.FUNCTIONAL_WEIGHTS)
        
    def InitializePhysicsSolutionStep(self):
        self._InitializeFluidSolutionStep()

    def _InitializeFluidSolutionStep(self):
        pass

    def InitializeAdjointPhysicsSolutionStep(self):
        self._InitializeAdjointFluidSolutionStep()

    def _InitializeAdjointFluidSolutionStep(self):
        self._InitializeAdjointFluidSolutionStepConvectionVelocity()
        self._InitializeAdjointPhysicsSolutionStepFunctionalDerivativeVelocity()

    def _InitializeAdjointFluidSolutionStepConvectionVelocity(self):
        convection_velocity = self.velocity_solutions[self.n_time_steps-self.time_step_counter,:,:].copy()
        self._GetAdjointSolver()._UpdateConvectionVelocityVariable(convection_velocity)
        self._GetMainModelPart().GetCommunicator().SynchronizeNonHistoricalVariable(KratosCFD.CONVECTION_VELOCITY)

    def _InitializeAdjointPhysicsSolutionStepFunctionalDerivativeVelocity(self):
        functional_derivative_velocity = self.velocity_solutions[self.n_time_steps-self.time_step_counter,:,:].copy()
        self._GetAdjointSolver()._UpdateFunctionalDerivativeVelocityVariable(functional_derivative_velocity)
        self._GetMainModelPart().GetCommunicator().SynchronizeNonHistoricalVariable(KratosMultiphysics.FUNCTIONAL_DERIVATIVE_VELOCITY)

    def FinalizeSolutionStep(self):
        super().FinalizeSolutionStep()
        self.EvaluateOptimizationRequiredGradients()
        if self.IsPhysicsStage():
            self.EvaluateFunctionalsInDeltaTime(self.time_step_counter-1)
        self.StoreTimeStepSolution()
    
    def StoreTimeStepSolution(self):
        if self.IsPhysicsStage():
            velocity = np.asarray(KratosMultiphysics.VariableUtils().GetSolutionStepValuesVector(self._GetLocalMeshNodes(), KratosMultiphysics.VELOCITY, 0, self.dim)).reshape(self.n_nodes, self.dim)
            pressure = np.asarray(KratosMultiphysics.VariableUtils().GetSolutionStepValuesVector(self._GetLocalMeshNodes(), KratosMultiphysics.PRESSURE, 0))
            self.velocity_solutions[self.time_step_counter-1,:,:] = velocity.copy()
            self.pressure_solutions[self.time_step_counter-1,:]   = pressure.copy()
        elif self.IsAdjointStage():
            velocity_adj = np.asarray(KratosMultiphysics.VariableUtils().GetSolutionStepValuesVector(self._GetLocalMeshNodes(), KratosMultiphysics.VELOCITY_ADJ, 0, self.dim)).reshape(self.n_nodes, self.dim)
            pressure_adj = np.asarray(KratosMultiphysics.VariableUtils().GetSolutionStepValuesVector(self._GetLocalMeshNodes(), KratosMultiphysics.PRESSURE_ADJ, 0))
            self.adjoint_velocity_solutions[self.n_time_steps-self.time_step_counter,:,:] = velocity_adj.copy()
            self.adjoint_pressure_solutions[self.n_time_steps-self.time_step_counter,:]   = pressure_adj.copy()

    def _ReInitializePhysics(self):
        """
        This method Re-Initializes the physical values in order to be consistent with the solution of the current optimization step
        """
        self.MpiPrint("--|" + self.topology_optimization_stage_str + "| RE-INITIALIZE PHYSICS")
        self.time = self.start_time
        self._GetSolver()._SetStartingTime(self.start_time)

    def _EvaluateOptimizationProblem(self, design_parameter = [], print_results = False):
        self.MpiPrint("\n--|EVALUATE OPTIMIZATION PROBLEM|")
        if (len(design_parameter) != 0):
            self._UpdateDesignParameterAndPhysicsParameters(design_parameter)
            self.EvaluateOptimizationRequiredGradients()
        self._EvaluateFunctionalAndDerivatives(print_results)
        self._EvaluateConstraintsAndDerivatives()
        if (print_results):
            self._PrintOptimizationProblem()
    
    def _UpdateOptimizationProblem(self, optimization_domain_design_parameter):
        design_parameter = self._InsertDesignParameterFromOptimizationDomain(optimization_domain_design_parameter)
        self._EvaluateOptimizationProblem(design_parameter, print_results=False)
        return self.functional, self._ExtractVariableInOptimizationDomain(self.functional_derivatives_wrt_design), self.constraints, self.constraints_derivatives_wrt_design
    
    def _EvaluateConstraintsAndDerivatives(self):
        self.MpiPrint("--|" + self.topology_optimization_stage_str + "| EVALUATE CONSTRAINTS AND DERIVATIVES")
        self._EvaluateVolumeConstraintAndDerivative()
        if (self.use_other_constraints):
            self._EvaluateOtherConstraintsAndDerivatives()
        self._UpdateVolumeConstraintDerivativesVariable()

    def _EvaluateOtherConstraintsAndDerivatives(self):
        if (self.use_wss_constraint):
            self._EvaluateWSSConstraintAndDerivative()

    def _EvaluateFunctionalAndDerivatives(self, print_functional=False):
        """
        This method is used to evaluate the functional value and derivatives w.r.t. the design parameter
        """
        self._SetTopologyOptimizationStage(3)
        self._EvaluateFunctional(print_functional)
        self._EvaluateFunctionalDerivatives()
    
    def _EvaluateResistanceFunctional(self, print_functional=False):
        """
        This method computes the resistance functional: int_{time}{int_{\Omega}{\\alpha||u||^2}}
        """
        self.fluid_functionals[0] = np.dot(self.resistance_functionals_in_delta_time, self.resistance_functional_time_steps_integration_weights)
        if self.IsMpiParallelism():
            self.fluid_functionals[0] = self.MpiSumLocalValues(self.fluid_functionals[0])
        if (self.first_iteration):
            self.initial_fluid_functionals_values[0] = self.fluid_functionals[0] 
        if (print_functional):
            self.MpiPrint("--|" + self.topology_optimization_stage_str + "| ---> Resistance Functional (no weight): " + str(self.fluid_functionals[0]))
        else:
            self.MpiPrint("--|" + self.topology_optimization_stage_str + "| ---> Resistance Functional")

    def _EvaluateResistanceFunctionalInDeltaTime(self):
        """
        This method computes the resistance functional for a single dt: int_{\Omega}{\\alpha||u||^2}
        """
        mp = self._GetPhysicsMainModelPart()
        velocity = np.asarray(KratosMultiphysics.VariableUtils().GetSolutionStepValuesVector(self._GetLocalMeshNodes(mp), KratosMultiphysics.VELOCITY, 0, self.dim)).reshape(self.n_nodes, self.dim)
        nodal_velocity_norm = np.linalg.norm(velocity, axis=1)
        integrand = self.resistance * (nodal_velocity_norm**2) #component-wise multiplication
        resistance_functional_in_delta_time = np.dot(self.nodal_domain_sizes, integrand)
        return resistance_functional_in_delta_time

    def _EvaluateStrainRateFunctional(self, print_functional=False):
        """
        This method computes the Strain-Rate functional: int_{time}{int_{\Omega}{\\2*mu*||1/2*[grad(u)+grad(u)^T]||^2}}
        """
        self.fluid_functionals[1] = np.dot(self.strain_rate_functionals_in_delta_time, self.strain_rate_functional_time_steps_integration_weights)
        if self.IsMpiParallelism():
            self.fluid_functionals[1] = self.MpiSumLocalValues(self.fluid_functionals[1])
        if (self.first_iteration):
            self.initial_fluid_functionals_values[1] = self.fluid_functionals[1]
        if (print_functional):
            self.MpiPrint("--|" + self.topology_optimization_stage_str + "| ---> Strain-Rate Functional (no weight): " + str(self.fluid_functionals[1]))
        else:
            self.MpiPrint("--|" + self.topology_optimization_stage_str + "| ---> Strain-Rate Functional")

    def _EvaluateStrainRateFunctionalInDeltaTime(self):
        """
        This method computes the resistance functional for a single dt: int_{\Omega}{\\2*mu*||1/2*[grad(u)+grad(u)^T]||^2}
        """
        vel_gradient_on_nodes = self._AssembleVelocityGradientOnNodes(buffer_id=0)
        vel_symmetric_gradient = 1.0/2.0 * (vel_gradient_on_nodes+(np.transpose(vel_gradient_on_nodes, axes=(0,2,1))))
        vel_symmetric_gradient_norm_squared = (np.linalg.norm(vel_symmetric_gradient, ord='fro', axis=(1, 2)))**2
        mu = self._GetViscosity()
        strain_rate_functional_in_delta_time = 2.0*mu* np.dot(vel_symmetric_gradient_norm_squared, self.nodal_domain_sizes)
        return strain_rate_functional_in_delta_time

    def _EvaluateVorticityFunctional(self, print_functional=False):
        """
        This method computes the Vorticity functional: int_{time}{int_{\Omega}{\\2*mu*||1/2*[grad(u)-grad(u)^T]||^2}}
        """
        self.fluid_functionals[2] = np.dot(self.vorticity_functionals_in_delta_time, self.vorticity_functional_time_steps_integration_weights)
        if self.IsMpiParallelism():
            self.fluid_functionals[2] = self.MpiSumLocalValues(self.fluid_functionals[2])
        if (self.first_iteration):
            self.initial_fluid_functionals_values[2] = self.fluid_functionals[2]
        if (print_functional):
            self.MpiPrint("--|" + self.topology_optimization_stage_str + "| ---> Vorticity Functional: (no weight) " + str(self.fluid_functionals[2])) 
        else:
            self.MpiPrint("--|" + self.topology_optimization_stage_str + "| ---> Vorticity Functional")

    def _EvaluateVorticityFunctionalInDeltaTime(self):
        """
        This method computes the resistance functional for a single dt: int_{\Omega}{\\2*mu*||1/2*[grad(u)-grad(u)^T]||^2}
        """
        vel_gradient_on_nodes = self._AssembleVelocityGradientOnNodes(buffer_id=0)
        vel_antisymmetric_gradient = 0.5 * (vel_gradient_on_nodes-(np.transpose(vel_gradient_on_nodes, axes=(0,2,1))))
        vel_antisymmetric_gradient_norm_squared = (np.linalg.norm(vel_antisymmetric_gradient, ord='fro', axis=(1, 2)))**2
        mu = self._GetViscosity()
        vorticity_functional_in_delta_time = 2.0*mu* np.dot(vel_antisymmetric_gradient_norm_squared, self.nodal_domain_sizes)
        return vorticity_functional_in_delta_time
    
    def _EvaluateFunctionalDerivatives(self):
        """
        This method is used to evaluate the functional derivatives w.r.t the design parameter
        """
        self._SetTopologyOptimizationStage(3)
        self.MpiPrint("--|" + self.topology_optimization_stage_str + "| EVALUATE FUNCTIONAL DERIVATIVES")
        self.functional_derivatives_wrt_design = np.asarray([self._ComputeFunctionalDerivatives()]).T
        self._UpdateFunctionalDerivativesVariable()
        if (self.first_iteration):
            self.initial_functional_derivatives_wrt_design = self.functional_derivatives_wrt_design

    def _ComputeFunctionalDerivatives(self):
        mask = self._GetOptimizationDomainNodesMask()
        temp_functional_derivatives_wrt_design = np.zeros(self.n_nodes)
        temp_functional_derivatives_wrt_design = self._ComputeFunctionalDerivativesFunctionalContribution()
        temp_functional_derivatives_wrt_design += self._ComputeFunctionalDerivativesPhysicsContribution()
        temp_functional_derivatives_wrt_design_projected = temp_functional_derivatives_wrt_design * self.design_parameter_projected_derivatives
        temp_functional_derivatives_wrt_design[mask] = self._ApplyDiffusiveFilterDerivative(temp_functional_derivatives_wrt_design_projected)
        return temp_functional_derivatives_wrt_design

    def _ComputeFunctionalDerivativesFunctionalContribution(self):
        return self._ComputeFunctionalDerivativesFluidFunctionalContribution()
    
    def _ComputeFunctionalDerivativesPhysicsContribution(self):
        return self._ComputeFunctionalDerivativesFluidPhysicsContribution()

    def _ComputeFunctionalDerivativesFluidFunctionalContribution(self):
        resistance_functional_derivatives = self._ComputeFunctionalDerivativesResistanceFunctionalContribution()
        fluid_functional_derivatives = self.functional_weights[0]*resistance_functional_derivatives
        return fluid_functional_derivatives 

    def _ComputeFunctionalDerivativesResistanceFunctionalContribution(self):
        resistance_functional_derivatives_in_delta_time = np.zeros((self.n_nodes, self.n_time_steps))
        for istep in range(self.n_time_steps):
            velocity = self.velocity_solutions[istep]
            resistance_functional_derivatives_in_delta_time[:,istep] = self.resistance_derivative_wrt_design_base * np.sum(velocity*velocity, axis=1) * self.nodal_domain_sizes 
        return np.einsum('ij,j->i', resistance_functional_derivatives_in_delta_time, self.resistance_functional_time_steps_integration_weights)

    def _ComputeFunctionalDerivativesFluidPhysicsContribution(self):
        fluid_physics_functional_derivatives_in_delta_time = np.zeros((self.n_nodes, self.n_time_steps))
        for istep in range(self.n_time_steps):
            velocity = self.velocity_solutions[istep]
            velocity_adjoint = self.adjoint_velocity_solutions[istep]       
            fluid_physics_functional_derivatives_in_delta_time[:,istep] = self.resistance_derivative_wrt_design_base * np.sum(velocity*velocity_adjoint, axis=1) * self.nodal_domain_sizes
        return np.einsum('ij,j->i', fluid_physics_functional_derivatives_in_delta_time, self.time_steps_integration_weights)

    def _EvaluateWSSConstraintAndDerivative(self):
        """
        For simplicity of notation we refer to the design parameter gradient as g.
        g: design parameter gradient
        gn: design parameter gradient norm 
        gn_int: integral over the domain of gn, used to normalize gn
        w: gn / gn_int (integral weights based on design parameter gradient)
        """
        g = np.asarray(KratosMultiphysics.VariableUtils().GetSolutionStepValuesVector(self._GetLocalMeshNodes(), KratosMultiphysics.DESIGN_PARAMETER_GRADIENT, 0, self.dim)).reshape(self.n_nodes, self.dim)
        gn = np.linalg.norm(g , axis=1)
        gn_max = np.max(gn)
        if (gn_max > 1e-14):
            gn_int = np.dot(gn, self.nodal_domain_sizes)
            v_grad = self._AssembleVelocityGradientOnNodes()
            nodal_tangents_to_g = self._EvaluateOrthogonalBasis(g)
            # prod = np.einsum('ijk,ik->ij',  nodal_tangents_to_g, g)
            psi_vect =np.einsum('ijk,ik,ilj->il', v_grad, g, nodal_tangents_to_g)
            psi =  np.linalg.norm(psi_vect, axis=1)
            mu = self._GetViscosity()
            self.wss_value = mu / gn_int * np.dot(psi, self.nodal_domain_sizes)
            self.wss_constraint = self.min_wss - self.wss_value
            self.constraints[self.wss_constraint_id] = self.wss_constraint
        else:
            self.wss_value = -1
            self.wss_constraint = -1
            self.constraints[self.wss_constraint_id] = self.wss_constraint
    
    ## UTILS
    def _SetTopologyOptimizationStage(self, problem_stage):
        """
        This method set the topology optimization problem phase:
        0: initialization, 1: NS solution, 2: ADJ NS solution, 3: Optimization Update
        """
        self.topology_optimization_stage = problem_stage 
        if self.CheckOptimizationStage(0):
            self.topology_optimization_stage_str = "INI"
        elif self.CheckOptimizationStage(1):
            self.topology_optimization_stage_str = "PHY"
        elif self.CheckOptimizationStage(2):
            self.topology_optimization_stage_str = "ADJ"
        elif self.CheckOptimizationStage(3):
            self.topology_optimization_stage_str = "OPT"
        else:
            self.topology_optimization_stage_str = "ERROR"
        self._GetSolver()._SetTopologyOptimizationStage(self.topology_optimization_stage)

    def GetTopologyOptimizationStage(self):
        """
        This method returns the topology optimization stage
        """
        return self.topology_optimization_stage

    def CheckOptimizationStage(self, check):
        """
        This method returns TRUE the topology optimization stage is 'check'
        """ 
        return (self.GetTopologyOptimizationStage() == check)
    
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
    
    def IsOptimizerStage(self):
        """
        This method returns TRUE the topology optimization stage is 3 <=> OPT
        """
        return self.CheckOptimizationStage(3) 
    
    def IsUnsteadySolution(self):
        return self._GetSolver()._IsUnsteady()
    
    def IssteadySolution(self):
        return not self.IsUnsteadySolution()
    
    def _GetComputingModelPart(self):
        """
        This method returns the computing model part for the current solver
        """    
        return self._GetSolver().GetComputingModelPart()
    
    def _GetMainModelPart(self):
        """
        This method returns the main model part for the current solver
        """    
        return self._GetSolver().GetMainModelPart()
    
    def _GetFluidModelPart(self):
        return self._GetMainModelPart()
    
    def _GetPhysicsMainModelPart(self):
        return self._GetPhysicsSolver().GetMainModelPart()
    
    def _GetAdjointMainModelPart(self):
        return self._GetAdjointSolver().GetMainModelPart()
    
    def _GetSubModelPart(self, super_model_part, sub_model_part_name):
        """
        This method returns the selected sub model part for the current solver
        """  
        if (super_model_part.HasSubModelPart(sub_model_part_name)):  
            return super_model_part.GetSubModelPart(sub_model_part_name)
        else:
            return None
    
    def _GetOptimizationDomain(self):
        """
        This method returns the optimization domain model part for the current solver
        """  
        opt_domain = self._GetSubModelPart(self._GetMainModelPart(), self.optimization_domain_name)
        if (opt_domain is None):
            return self._GetMainModelPart()
        else:
            return opt_domain
    
    def _GetNonOptimizationDomain(self):
        """
        This method returns the non-optimization domain model part for the current solver
        """  
        return self._GetSubModelPart(self._GetMainModelPart(), self.non_optimization_domain_name)

    def _GetNodes(self):
        """
        This method returns the nodes for the current solver
        """
        return self._GetComputingModelPart().Nodes
    
    def _GetElements(self):
        """
        This method returns the elements for the current solver
        """
        return self._GetComputingModelPart().Elements
    
    def _GetSimulationName(self):
        """
        This method returns simulation name
        """
        return "--|" + self.topology_optimization_stage_str + "| Solution Analysis"
    
    def CustomProblem(self):
        """
        This method is used to do custom actions in problem tests
        """
        dumb = 0
        return dumb
    
    def _GetNodeCoordinates(self, node):
        if (self.dim == 2):
            return np.asarray([node.X, node.Y])
        else:
            return np.asarray([node.X, node.Y, node.Z])

    def _GetNodesSetCoordinates(self, nodes):
        points = np.zeros((len(nodes),self.dim))
        count = 0
        for node in nodes:
            points[count,:] = self._GetNodeCoordinates(node)
            count +=1
        return points

    def _GetOptimizationDomainNodesMask(self):
        return self.optimization_domain_nodes_mask
    
    ## PRINTS
    def _PrintOptimizationProblem(self):
        self.MpiPrint("\n--|PRINT OPTIMIZATION PROBLEM DATA|", min_echo=0)
        self._PrintFunctionals()
        self._PrintConstraints()
    
    def _PrintFunctionals(self):
        self._PrintTotalFunctional()
        self._PrintFluidFunctionals()

    def _PrintTotalFunctional(self):
        self.MpiPrint("--|" + self.topology_optimization_stage_str + "| TOTAL FUNCTIONAL  : " +  str(self.functional), min_echo=0)
        self.MpiPrint("--|" + self.topology_optimization_stage_str + "| INITIAL FUNCTIONAL: " +  str(self.initial_functional), min_echo=0)
        
    def _PrintFluidFunctionals(self):
        if (abs(self.functional_weights[0]) > 1e-10):
            self.MpiPrint("--|" + self.topology_optimization_stage_str + "| ---> Resistance Functional (" + str(self.functional_weights[0]) + "): " + str(self.weighted_functionals[0]), min_echo=0)
        if (abs(self.functional_weights[1]) > 1e-10):
            self.MpiPrint("--|" + self.topology_optimization_stage_str + "| ---> Strain-Rate Functional (" + str(self.functional_weights[1]) + "): " + str(self.weighted_functionals[1]), min_echo=0)
        if (abs(self.functional_weights[2]) > 1e-10):
            self.MpiPrint("--|" + self.topology_optimization_stage_str + "| ---> Vorticity Functional (" + str(self.functional_weights[2]) + "): " + str(self.weighted_functionals[2]), min_echo=0)
    
    def _PrintFunctionalsToFile(self):
        if self.MpiRunOnlyRank(0):
                destination_path = self.optimization_parameters["optimization_settings"]["solution_output_settings"]["output_path"].GetString()
                destination_path += "/files"
                with open(destination_path+"/functional_history.txt", "a") as file:
                    file.write(str(self.functional) + " ")
                    for ifunc in range(self.n_functionals): 
                        file.write(str(self.weighted_functionals[ifunc]) + " ")
                    file.write("\n")

    def _EvaluateDesignParameterChange(self):
            design_parameter_converged = False
            design_parameter_per_rank = self.data_communicator.AllGathervDoubles(self.design_parameter)
            old_design_parameter_per_rank = self.data_communicator.AllGathervDoubles(self.old_design_parameter)
            if (self.MpiRunOnlyRank(0)):
                design_parameter = np.concatenate(design_parameter_per_rank)
                old_design_parameter = np.concatenate(old_design_parameter_per_rank)
                old_design_parameter_norm = np.linalg.norm(1-old_design_parameter)
                self.design_parameter_change = np.linalg.norm(design_parameter-old_design_parameter)
                if (old_design_parameter_norm > 1e-10):
                    # change evaluated on the amount of fluid, that's why we have ||(1-design)-(1-old_design)|| / ||1-old_design||
                    self.design_parameter_change /= old_design_parameter_norm
                else:
                    # change evaluated on the amount of fluid, that's why we have ||design-old_design|| / ||old_design||
                    old_design_parameter_norm = np.linalg.norm(self.old_design_parameter)
                    self.design_parameter_change /= old_design_parameter_norm
                if (self.opt_it > 1):
                    design_parameter_converged = (self.design_parameter_change < self.design_parameter_change_toll)
                else:
                    design_parameter_converged = False
            design_parameter_converged = self.data_communicator.Broadcast(design_parameter_converged, source_rank=0)
            self.design_parameter_change = self.data_communicator.Broadcast(self.design_parameter_change, source_rank=0)
            self.MpiPrint("--|" + self.topology_optimization_stage_str + "| DESIGN PARAMETER CHANGE: " + str(self.design_parameter_change))
            return design_parameter_converged       

    def _ResetFunctionalOutput(self):
        if self.MpiRunOnlyRank(0):
            dest_dir = Path(self.optimization_parameters["optimization_settings"]["solution_output_settings"]["output_path"].GetString()) / "files"
            dest_dir.mkdir(parents=True, exist_ok=True)
            destination_path = self.optimization_parameters["optimization_settings"]["solution_output_settings"]["output_path"].GetString()
            destination_path += "/files"
            with open(destination_path+"/functional_history.txt", "w") as file:
                file.write("TOTAL | RESISTANCE | STRAIN_RATE | 2*VORTICITY | OUTLET_CONCENTRATION | REGION_TEMPERATURE | DIFFUSIVE_TRANSPORT | CONVECTIVE_TRASNPORT | REACTIVE_TRASNPORT | SOURCE_TRANSPORT\n")
                file.write("1 ")
                for ifunc in range(self.n_functionals):
                    file.write(str(self.functional_weights[ifunc]) + " ")
                file.write("\n")

    def _PrintConstraints(self):
        self._PrintVolumeConstraint()
        if (self.use_other_constraints):
            self._PrintOtherConstraints()

    def _PrintOtherConstraints(self):
        if (self.use_wss_constraint):
            self._PrintWSSConstraint()
    
    def _PrintVolumeConstraint(self):
        self.MpiPrint("--|" + self.topology_optimization_stage_str + "| VOLUME FRACTION: " + str(self.volume_fraction), min_echo=0)
        self.MpiPrint("--|" + self.topology_optimization_stage_str + "| ---> Volume Constraint: " + str(self.volume_constraint), min_echo=0)

    def _PrintWSSConstraint(self):
        self.MpiPrint("--|" + self.topology_optimization_stage_str + "| WSS VALUE: " + str(self.wss_value), min_echo=0)
        # self.MpiPrint("--|" + self.topology_optimization_stage_str + "| ---> WSS Resistance: " + str(self.resistance_parameters["value_full"].GetDouble()), min_echo=0)
        self.MpiPrint("--|" + self.topology_optimization_stage_str + "| ---> WSS Constraint: " + str(self.wss_constraint), min_echo=0)

    def PrintAnalysisStageProgressInformation(self):
        KratosMultiphysics.Logger.PrintInfo(self._GetSimulationName(), "TOTAL STEP: ", self._GetComputingModelPart().ProcessInfo[KratosMultiphysics.STEP])
        if self.IsPhysicsStage(): # NS
            KratosMultiphysics.Logger.PrintInfo(self._GetSimulationName(), "NS STEP: ", self._GetComputingModelPart().ProcessInfo[KratosCFD.FLUID_TOP_OPT_NS_STEP])
        elif self.IsAdjointStage(): # ADJ
            KratosMultiphysics.Logger.PrintInfo(self._GetSimulationName(), "ADJ NS STEP: ", self._GetComputingModelPart().ProcessInfo[KratosCFD.FLUID_TOP_OPT_ADJ_NS_STEP])
        else:
            KratosMultiphysics.Logger.PrintInfo(self._GetSimulationName(), "Invalid Value of the Topology Optimization Stage. TOP_OPT_STAGE: ", self.topology_optimization_stage)
        KratosMultiphysics.Logger.PrintInfo(self._GetSimulationName(), "TIME: ", self.time)

    def _ComputeDesignParameterDiffusiveFilterUtilities(self):
        self._InitializeDiffusiveFilter()
        self.design_parameter_filtered = np.zeros(self.n_nodes)

    def _InitializeDiffusiveFilter(self):
        self.diffusive_filter_settings = self.optimization_settings["diffusive_filter_settings"]
        self.apply_diffusive_filter = self.diffusive_filter_settings["use_filter"].GetBool()
        self.diffusive_filter_type = self.diffusive_filter_settings["filter_type"].GetString()
        self.diffusive_filter_radius = self.diffusive_filter_settings["radius"].GetDouble()
        self.diffusive_filter_type_settings = self.diffusive_filter_settings["type_settings"]
        if self.apply_diffusive_filter:
            if self.diffusive_filter_type == "pde":
                self._InitializePdeDiffusiveFilter()
            else:
                self._InitializeDiscreteDiffusiveFilter()

    def _ComputeDesignParameterProjectiveFilterUtilities(self):
        self._InitializeProjectiveFilter()
        self.design_parameter_projected = np.zeros(self.n_nodes)
        self.design_parameter_projected_derivatives = np.ones(self.n_nodes)

    def _InitializeProjectiveFilter(self):
        self.projective_filter_settings = self.optimization_settings["projective_filter_settings"]
        self.apply_projective_filter = self.projective_filter_settings["use_filter"].GetBool()
        self.projective_filter_min_mean = self.projective_filter_settings["min_max_mean"][0].GetDouble()
        self.projective_filter_max_mean = self.projective_filter_settings["min_max_mean"][1].GetDouble()
        self.projective_filter_min_projection_slope = self.projective_filter_settings["min_max_projection_slope"][0].GetDouble()
        self.projective_filter_max_projection_slope = self.projective_filter_settings["min_max_projection_slope"][1].GetDouble()
        self.projective_filter_activation_change = self.projective_filter_settings["activation_change"].GetDouble()
        self.projective_filter_slope = self.projective_filter_min_projection_slope

    def _InitializeDiscreteDiffusiveFilter(self):
        if (self.diffusive_filter_radius < 1e-10): #ensures that if no filter is imposed, at least the node itself is in neighboring nodes
            self.diffusive_filter_radius = 1e-10
        self.MpiPrint("--|" + self.topology_optimization_stage_str + "| ---> Initialize Discrete Diffusive Filter")
        mp = self._GetComputingModelPart()
        mask = self._GetOptimizationDomainNodesMask()
        global_ids = [self.nodes_ids_local_partition_to_global_dictionary[int(i)] for i in mask]
        only_opt_mp_nodes = self._GetModelPartNodesSubset(mp, global_ids)
        points = self._GetNodesSetCoordinates(only_opt_mp_nodes)
        nodes_tree = KDTree(points)
        self.nodes_connectivity_matrix = nodes_tree.sparse_distance_matrix(nodes_tree, self.diffusive_filter_radius, output_type="dok_matrix").tocsr()
        self.nodes_connectivity_matrix.data *= -1
        self.nodes_connectivity_matrix.data += self.diffusive_filter_radius
        self._CorrectNodesConnectivityMatrixWithSymmetry()
        self.nodes_connectivity_weights_sum = np.array(self.nodes_connectivity_matrix.sum(axis=1)).flatten()  # Sum of each row as a 1D array
        # Normalization step: Divide non-zero entries by the corresponding row sum
        row_indices = np.repeat(np.arange(self.nodes_connectivity_matrix.shape[0]), np.diff(self.nodes_connectivity_matrix.indptr))
        self.nodes_connectivity_matrix.data /= self.nodes_connectivity_weights_sum[row_indices]
        self.nodes_connectivity_matrix_for_derivatives = self.nodes_connectivity_matrix.copy()
        self._CorrectNodesConnectivityMatrixForDerivativesWithSymmetryAndTranspose()

    def _InitializePdeDiffusiveFilter(self):
        self.MpiPrint("--|" + self.topology_optimization_stage_str + "| ---> Initialize PDE Diffusive Filter")
        self.pde_diffusive_filter_process = ApplyTopologyOptimizationPdeFilterProcess(self.model, self.diffusive_filter_type_settings, self.diffusive_filter_radius, self._GetMainModelPart(),  self._GetOptimizationDomain(), self._GetOptimizationDomainNodesMask(), self.nodes_ids_global_to_local_partition_dictionary)

    def _CorrectNodesConnectivityMatrixWithSymmetry(self):
        """
        Correct the nodes connnectvity matrix to take into account the node that would be in the other side of the symmetry axis
        Note: IT WORKS WELL WITH REGULAR MESHES AND FILTER RADIUS OF THE SIZE OF THE ELEMENTS
        """
        if (self.symmetry_enabled):
            symm_to_opt_mask_nodes = self.global_to_opt_mask_nodes[self.symmetry_nodes_mask]
            for inode in range(symm_to_opt_mask_nodes.size):
                mask_node = symm_to_opt_mask_nodes[inode]
                if (mask_node != -1):
                    row = mask_node
                    start = self.nodes_connectivity_matrix.indptr[row]
                    end = self.nodes_connectivity_matrix.indptr[row + 1]
                    for idx in range(start, end):
                        col = self.nodes_connectivity_matrix.indices[idx]
                        if col != row:
                            self.nodes_connectivity_matrix.data[idx] *= 2.0

    def _CorrectNodesConnectivityMatrixForDerivativesWithSymmetryAndTranspose(self):
        """
        Correct the nodes connnectvity matrix for drivatives to not take intop account a oduble value of the gradient in a node
        Note: IT WORKS WELL WITH REGULAR MESHES AND FILTER RADIUS OF THE SIZE OF THE ELEMENTS
        """
        if (self.symmetry_enabled):
            symm_to_opt_mask_nodes = self.global_to_opt_mask_nodes[self.symmetry_nodes_mask]
            for inode in range(symm_to_opt_mask_nodes.size):
                mask_node = symm_to_opt_mask_nodes[inode]
                if (mask_node != -1):
                    row = mask_node
                    start = self.nodes_connectivity_matrix_for_derivatives.indptr[row]
                    end = self.nodes_connectivity_matrix_for_derivatives.indptr[row + 1]
                    for idx in range(start, end):
                        col = self.nodes_connectivity_matrix_for_derivatives.indices[idx]
                        if col != row:
                            self.nodes_connectivity_matrix_for_derivatives.data[idx] /= 2.0
        self.nodes_connectivity_matrix_for_derivatives = self.nodes_connectivity_matrix_for_derivatives.transpose().tocsr()
        if (self.symmetry_enabled):
            symm_to_opt_mask_nodes = self.global_to_opt_mask_nodes[self.symmetry_nodes_mask]
            for inode in range(symm_to_opt_mask_nodes.size):
                mask_node = symm_to_opt_mask_nodes[inode]
                if (mask_node != -1):
                    row = mask_node
                    start = self.nodes_connectivity_matrix_for_derivatives.indptr[row]
                    end = self.nodes_connectivity_matrix_for_derivatives.indptr[row + 1]
                    for idx in range(start, end):
                        col = self.nodes_connectivity_matrix_for_derivatives.indices[idx]
                        if col != row:
                            self.nodes_connectivity_matrix_for_derivatives.data[idx] *= 2.0
        
    def _ApplyDesignParameterDiffusiveFilter(self, design_parameter):
        self.MpiPrint("--|" + self.topology_optimization_stage_str + "| --> Apply Diffusive Filter: " + str(self.apply_diffusive_filter))
        self.design_parameter_filtered = design_parameter
        mask = self._GetOptimizationDomainNodesMask()
        self.design_parameter_filtered[mask] = self._ApplyDiffusiveFilter(design_parameter)

    def _ApplyDesignParameterProjectiveFilter(self, design_parameter):
        self.MpiPrint("--|" + self.topology_optimization_stage_str + "| --> Apply Projective Filter: " + str(self.apply_projective_filter))
        mask = self._GetOptimizationDomainNodesMask()
        design_parameter = np.clip(design_parameter, a_min=0.0, a_max=1.0)
        self.design_parameter_projected = design_parameter
        self.design_parameter_projected_derivatives = np.ones(self.n_nodes)
        if ((self.opt_it > 10) and (self.apply_projective_filter)):
            design_change_ratio = 1.0-min(1.0, max(0.0, self.design_parameter_change-self.design_parameter_change_toll)/(self.projective_filter_activation_change-self.design_parameter_change_toll))
            if (abs(design_change_ratio) < 1e-15): # design_change_ratio == 0
                new_projection_slope = self.projective_filter_min_projection_slope
            else:
                self.MpiPrint("--|" + self.topology_optimization_stage_str + "| --> Apply Projective Filter: ---> APPLY PROJECTION")
                if ((abs(1-design_change_ratio) < 1e-15)): # design_change_ratio == 1
                    new_projection_slope = self.projective_filter_max_projection_slope
                else: # design_change_ratio \in (0,1)
                    new_projection_slope = self.projective_filter_min_projection_slope + (self.projective_filter_max_projection_slope-self.projective_filter_min_projection_slope)*design_change_ratio    
            if (new_projection_slope > self.projective_filter_slope):
                self.projective_filter_slope = new_projection_slope
            design_parameter_mean = min(max(self.projective_filter_min_mean,1.0-self.volume_fraction),self.projective_filter_max_mean)
            factor = np.tanh(self.projective_filter_slope*design_parameter_mean) + np.tanh(self.projective_filter_slope*(1-design_parameter_mean))
            constant = np.tanh(self.projective_filter_slope*design_parameter_mean)
            projection_argument = self.projective_filter_slope*(design_parameter-design_parameter_mean)
            self.design_parameter_projected[mask] = ((constant + np.tanh(projection_argument)) / factor)[mask]
            self.design_parameter_projected_derivatives[mask] = ((self.projective_filter_slope / factor) / (np.cosh(projection_argument)**2))[mask]
        else:
            self.design_parameter_projected = np.clip(self.design_parameter_projected, a_min=0.0, a_max=1.0)
                
    def _ApplyDiffusiveFilter(self, scalar_variable):
        scalar_variable_in_opt_domain = scalar_variable[self._GetOptimizationDomainNodesMask()]
        if self.apply_diffusive_filter:
            if (self.diffusive_filter_type == "pde"):
                self.MpiPrint("--|" + self.topology_optimization_stage_str + "| ----> PDE Filter for Design Parameter")
                self._InitializePdeDiffusiveFilterExecution(scalar_variable)
                self.pde_diffusive_filter_process.Execute()
                return self._FinalizePdeDiffusiveFilterExecution()
            else:
                self.MpiPrint("--|" + self.topology_optimization_stage_str + "| ----> Discrete Filter for Design Parameter")
                return (self.nodes_connectivity_matrix @ scalar_variable_in_opt_domain)
        else:
            return scalar_variable_in_opt_domain
        
    def _InitializePdeDiffusiveFilterExecution(self, pde_filter_source):
        self.pde_diffusive_filter_process._UpdatePdeFilterSourceTerm(pde_filter_source)

    def _FinalizePdeDiffusiveFilterExecution(self):
        return self.pde_diffusive_filter_process._GetLastFilteredValue()
    
    def _ApplyDiffusiveFilterDerivative(self, scalar_variable_derivative):
        scalar_variable_derivative_in_opt_domain = scalar_variable_derivative[self._GetOptimizationDomainNodesMask()]
        if self.apply_diffusive_filter:
            if (self.diffusive_filter_type == "pde"):
                self.MpiPrint("--|" + self.topology_optimization_stage_str + "| ----> PDE Filter for Functional Derivative")
                self._InitializePdeDiffusiveFilterExecution(scalar_variable_derivative)
                self.pde_diffusive_filter_process.Execute()
                return self._FinalizePdeDiffusiveFilterExecution()
            else:
                self.MpiPrint("--|" + self.topology_optimization_stage_str + "| ----> Discrete Filter for Functional Derivative")
                return self.nodes_connectivity_matrix_for_derivatives @ scalar_variable_derivative_in_opt_domain
        else:
            return scalar_variable_derivative_in_opt_domain
        

    def _UpdateOptimizationProblemPhysics(self, optimization_domain_design_parameter):
        self.design_parameter = self._InsertDesignParameterFromOptimizationDomain(optimization_domain_design_parameter)
        self._SolveTopologyOptimizationStepPhysics()
        self._EvaluateOptimizationProblem(print_results=True)
        return self.functional, self._ExtractVariableInOptimizationDomain(self.functional_derivatives_wrt_design), self.constraints, self.constraints_derivatives_wrt_design
        
    def _GetDesignParameterVariable(self):
        design_parameter = np.zeros(self.n_nodes)
        count = 0
        for node in self._GetLocalMeshNodes():
            design_parameter[count] = node.GetSolutionStepValue(KratosMultiphysics.DESIGN_PARAMETER)
            count += 1
        return design_parameter
    
    def OptimizationOutputSolutionStep(self):
        self.MpiPrint("\n--| PRINT OPTIMIZATION STEP OUTPUT TO FILES")
        """This function printed / writes output files after the solution of a optimization step
        """
        execute_was_called = False
        for output_process in self._GetListOfOptimizationOutputProcesses():
            if output_process.IsOutputStep():
                if not execute_was_called:
                    for process in self._GetListOfProcesses():
                        process.ExecuteBeforeOutputStep()
                    execute_was_called = True
                output_process.PrintOutput()
        if execute_was_called:
            for process in self._GetListOfProcesses():
                process.ExecuteAfterOutputStep()

    def TimeOutputSolutionStep(self, print_time_step=True):
        """This function printed / writes output files after the solution of a time step
        """
        if (print_time_step):
            self.MpiPrint("\n--| PRINT TIME SOLUTION OUTPUT TO FILES")
            for istep in range(self.n_time_steps):
                self.SetTimeOutputSolutionStepPhysicsVariables(istep)
                for output_process in self._GetListOfTimeOutputProcesses():
                    self.PrintTimeOutputProcess(output_process, istep+1)

    def SetTimeOutputSolutionStepPhysicsVariables(self, time_step_id):
        self._SetTimeOutputSolutionStepFluidVariables(time_step_id)

    def _SetTimeOutputSolutionStepFluidVariables(self, time_step_id): 
        KratosMultiphysics.VariableUtils().SetSolutionStepValuesVector(self._GetLocalMeshNodes(), KratosMultiphysics.VELOCITY, self.velocity_solutions[time_step_id].flatten(), 0)
        KratosMultiphysics.VariableUtils().SetSolutionStepValuesVector(self._GetLocalMeshNodes(), KratosMultiphysics.PRESSURE, self.pressure_solutions[time_step_id], 0)  
        KratosMultiphysics.VariableUtils().SetSolutionStepValuesVector(self._GetLocalMeshNodes(), KratosMultiphysics.VELOCITY_ADJ, self.adjoint_velocity_solutions[time_step_id].flatten(), 0)  
        KratosMultiphysics.VariableUtils().SetSolutionStepValuesVector(self._GetLocalMeshNodes(), KratosMultiphysics.PRESSURE_ADJ, self.adjoint_pressure_solutions[time_step_id], 0) 

    def PrintTimeOutputProcess(self, output_process, time_step_id):
        time_step_counter = str(time_step_id).zfill(len(str(self.n_time_steps))) 
        file_name = "time_step_" + time_step_counter
        if type(output_process).__name__ == "VtkOutputProcess":
            output_process.vtk_io.PrintOutput(file_name)
        elif type(output_process).__name__ == "VtuOutputProcess":
            for vtu_output in output_process.vtu_output_ios:
                vtu_output.PrintOutput(f"{output_process.output_path / file_name}")

    def _GetListOfOptimizationOutputProcesses(self):
        """This function returns the list of output processes involved in this Analysis
        """
        if not hasattr(self, '_list_of_optimization_output_processes'):
            raise Exception("The list of optimization output-processes was not yet created!")
        return self._list_of_optimization_output_processes
    
    def _GetListOfTimeOutputProcesses(self):
        """This function returns the list of output processes involved in this Analysis
        """
        if not hasattr(self, '_list_of_optimization_output_processes'):
            raise Exception("The list of optimization output-processes was not yet created!")
        return self._list_of_time_output_processes

    def _PrintOptimizationSolution(self):
        self.OptimizationOutputSolutionStep()
        self._CorrectPvtuFilesInOptimizationVtuOutput()
        if self._IsTopologyOptimizationSolutionEnd():
            self.TimeOutputSolutionStep()
            self._CorrectPvtuFilesInTimeVtuOutput()
        self._PrintFunctionalsToFile()
        if self.first_iteration:
            self._CopyInputFiles()
        
    def _CopyInputFiles(self):
        if self.MpiRunOnlyRank(0):
            destination_path = self.optimization_parameters["optimization_settings"]["solution_output_settings"]["output_path"].GetString()
            destination_path += "/files"
            self._CopyFile(src_path="ProjectParameters.json", dst_path=destination_path+"/ProjectParameters.json")
            self._CopyFile(src_path="OptimizationParameters.json", dst_path=destination_path+"/OptimizationParameters.json")

    def _CopyFile(self, src_path, dst_path):
        src = Path(src_path)
        dst = Path(dst_path)
        if not src.exists():
            raise FileNotFoundError(f"Source file not found: {src}")
        dst.parent.mkdir(parents=True, exist_ok=True)  # ensure destination folder exists
        shutil.copy2(src, dst)

    def _InitializeRemeshing(self):
        self.remeshing_settings = self.optimization_settings["remeshing_settings"]
        self.enable_remeshing = self.remeshing_settings["remesh"].GetBool()
        self.remeshing_levelset = self.remeshing_settings["level_set"].GetDouble()
        self.remeshing_min_element_size = self.remeshing_settings["min_max_element_size"][0].GetDouble()
        self.remeshing_max_element_size = self.remeshing_settings["min_max_element_size"][1].GetDouble()
        if (self.IsRemeshingEnabled()):
            self.MpiPrint("--|" + self.topology_optimization_stage_str + "| INITIALIZE REMESHING PROCESS")
            main_mp = self._GetMainModelPart()
            # Create find_nodekl_h process
            self.find_nodal_h = KratosMultiphysics.FindNodalHNonHistoricalProcess(main_mp)
            self.find_nodal_h.Execute()

            # Create Local Hessian and Remeshing processes
            metric_param = KratosMultiphysics.Parameters(
                                                """{
                                                    "hessian_strategy_parameters"              :{
                                                            "estimate_interpolation_error"     : false,
                                                            "interpolation_error"              : 0.004,
                                                            "mesh_dependent_constant"          : 0.1
                                                    },
                                                    "minimal_size"                      : """ + str(self.remeshing_min_element_size) + """,
                                                    "maximal_size"                      : """ + str(self.remeshing_max_element_size) + """,
                                                    "enforce_current"                   : false,
                                                    "anisotropy_remeshing"              : true,
                                                    "enforced_anisotropy_parameters":{
                                                        "reference_variable_name"          : "DISTANCE",
                                                        "hmin_over_hmax_anisotropic_ratio" : 0.15,
                                                        "boundary_layer_max_distance"      : 1.0,
                                                        "interpolation"                    : "Linear"
                                                    }
                                                }"""
                                            )
            remesh_param = KratosMultiphysics.Parameters("""{ }""")
            if (self.dim == 2):
                self.local_hessian = KratosMultiphysics.ComputeHessianSolMetricProcess2D(main_mp, KratosMultiphysics.DISTANCE, metric_param)
                self.mmg_process = KratosMMG.MmgProcess2D(main_mp, remesh_param)
            elif (self.dim == 3):
                self.local_hessian = KratosMMG.ComputeHessianSolMetricProcess3D(main_mp, KratosMultiphysics.DISTANCE, metric_param)
                self.mmg_process = KratosMMG.MmgProcess3D(main_mp, remesh_param)
            self.local_hessian.Execute()

    def IsRemeshingEnabled(self):
        return self.enable_remeshing

    def _Remesh(self):
        if (self.enable_remeshing):
            self.MpiPrint("--|" + self.topology_optimization_stage_str + "| DOMAIN REMESHING")
            self.find_nodal_h.Execute()
            self.local_hessian.Execute()
            self.mmg_process.Execute()
            self._PreprocessRemeshedGeometry()

    def _PreprocessRemeshedGeometry(self):
        # GEOMETRICAL PREPROCESSING
        self._CreateNodesIdsDictionary()
        self._CreateElementsIdsDictionary()
        self._ComputeDomainSize()

        # OPTIMIZATION PREPROCESSING
        self._ComputeDesignParameterFilterUtilities()
        self._ResetConstraints()
        self.design_parameter = self._GetDesignParameterVariable()
        self._ResetPhysicsParameters()
        self._UpdatePhysicsParameters()
        self._PreprocessDerivatives() 

    def _GetPhysicsMainModelPartsList(self):
        return [self._GetPhysicsSolver().main_model_part]
    
    def EvaluateFunctionals(self, print_functional):
        self.EvaluateOptimizationRequiredGradients()
        self.EvaluatePhysicsFunctionals(print_functional)        

    def EvaluatePhysicsFunctionals(self, print_functional):
        self.EvaluateFluidFunctionals(print_functional)

    def EvaluateFluidFunctionals(self, print_functional):
        if (abs(self.normalized_fluid_functional_weights[0]) > 1e-10):
            self._EvaluateResistanceFunctional(print_functional)
        if (abs(self.normalized_fluid_functional_weights[1]) > 1e-10):
            self._EvaluateStrainRateFunctional(print_functional)
        if (abs(self.normalized_fluid_functional_weights[2]) > 1e-10):
            self._EvaluateVorticityFunctional(print_functional)

    def EvaluateFunctionalsInDeltaTime(self, time_step_id):
        self.EvaluatePhysicsFunctionalInDeltaTime(time_step_id)

    def EvaluatePhysicsFunctionalInDeltaTime(self, time_step_id):
        self.EvaluateFluidFunctionalsInDeltaTime(time_step_id)

    def EvaluateFluidFunctionalsInDeltaTime(self, time_step_id):
        if (abs(self.normalized_fluid_functional_weights[0]) > 1e-10):
            self.resistance_functionals_in_delta_time[time_step_id]  = self._EvaluateResistanceFunctionalInDeltaTime()
        if (abs(self.normalized_fluid_functional_weights[1]) > 1e-10):
            self.strain_rate_functionals_in_delta_time[time_step_id] = self._EvaluateStrainRateFunctionalInDeltaTime()
        if (abs(self.normalized_fluid_functional_weights[2]) > 1e-10):
            self.vorticity_functionals_in_delta_time[time_step_id]   = self._EvaluateVorticityFunctionalInDeltaTime()

    def EvaluateTotalFunctional(self):
        self.functionals = np.concatenate((self.fluid_functionals, np.zeros(self.n_functionals-self.n_fluid_functionals)))
        self.weighted_functionals = self.functional_weights * self.functionals
        self.functional = np.sum(self.weighted_functionals)
        if (self.first_iteration):
            self.initial_functional = self.functional

    def _CheckMaterialProperties(self, check = False):
        if (check):
            self.MpiPrint("--|CHECK| Check Fluid Properties")
            self._GetSolver()._CheckMaterialProperties()

    def _UpdateRelevantPhysicsVariables(self):
        pass

    def _UpdateRelevantAdjointVariables(self):
        pass

    def EvaluateOptimizationRequiredGradients(self):
        # evaluate gradients used for functional and constraints evaluation and their derivatives
        self._EvaluateVelocityGradient()
        # WSS constraint
        self._ComputeScalarVariableNodalGradient(KratosMultiphysics.DESIGN_PARAMETER, KratosMultiphysics.DESIGN_PARAMETER_GRADIENT)

    def _EvaluateVelocityGradient(self):
        self._ComputeScalarVariableNodalGradient(KratosMultiphysics.VELOCITY_X, KratosMultiphysics.VELOCITY_X_GRADIENT)
        self._ComputeScalarVariableNodalGradient(KratosMultiphysics.VELOCITY_Y, KratosMultiphysics.VELOCITY_Y_GRADIENT)
        if (self.dim == 3):
            self._ComputeScalarVariableNodalGradient(KratosMultiphysics.VELOCITY_Z, KratosMultiphysics.VELOCITY_Z_GRADIENT)

    def _AssembleVelocityGradientOnNodes(self, buffer_id):
        gradient_x = np.asarray(KratosMultiphysics.VariableUtils().GetSolutionStepValuesVector(self._GetLocalMeshNodes(), KratosMultiphysics.VELOCITY_X_GRADIENT, buffer_id, self.dim)).reshape(self.n_nodes, self.dim)
        gradient_y = np.asarray(KratosMultiphysics.VariableUtils().GetSolutionStepValuesVector(self._GetLocalMeshNodes(), KratosMultiphysics.VELOCITY_Y_GRADIENT, buffer_id, self.dim)).reshape(self.n_nodes, self.dim)
        gradient = np.zeros((self.n_nodes, self.dim, self.dim))
        gradient[:,0,:] = gradient_x
        gradient[:,1,:] = gradient_y
        if (self.dim == 3):
            gradient_z = np.asarray(KratosMultiphysics.VariableUtils().GetSolutionStepValuesVector(self._GetLocalMeshNodes(), KratosMultiphysics.VELOCITY_Z_GRADIENT, buffer_id, self.dim)).reshape(self.n_nodes, self.dim)
            gradient[:,2,:] = gradient_z
        return gradient
    
    def _EvaluateOrthogonalBasis(self, vector_in_nodes):
        mask = np.linalg.norm(vector_in_nodes, axis=1) > 1e-10 # relevan nodes mask. Avoids division by zero
        orthogonal_basis = np.zeros((self.n_nodes,self.dim-1,self.dim))
        orthogonal_basis[mask,0,0] = -vector_in_nodes[mask,1]
        orthogonal_basis[mask,0,1] =  vector_in_nodes[mask,0]
        orthogonal_basis_0_norm = np.linalg.norm(orthogonal_basis[:,0,:], axis=1)
        orthogonal_basis[mask,0,:] /= orthogonal_basis_0_norm[mask, np.newaxis]
        if (self.dim == 3):
            orthogonal_basis[mask,1,0] =  -vector_in_nodes[mask,2] + (vector_in_nodes[mask,1]*vector_in_nodes[mask,1]*vector_in_nodes[mask,2]) / (orthogonal_basis_0_norm[mask]**2)
            orthogonal_basis[mask,1,1] = -(vector_in_nodes[mask,0]*vector_in_nodes[mask,1]*vector_in_nodes[mask,2]) / orthogonal_basis_0_norm[mask]**2
            orthogonal_basis[mask,1,2] =   vector_in_nodes[mask,0]
            orthogonal_basis_1_norm = np.linalg.norm(orthogonal_basis[:,1,:], axis=1)
            orthogonal_basis[mask,1,:] /= orthogonal_basis_1_norm[mask, np.newaxis]
        return orthogonal_basis

    def GetDefaultPhysicsParametersSettings(self):
        ##settings string in json format
        default_physics_parameters_settings = KratosMultiphysics.Parameters("""
        {
            "resistance": [{
                "interpolation_method": "hyperbolic",
                "value_void"        : 0.0,
                "value_full"        : 1e4,
                "interpolation_slope": {
                    "initial_slope": 1.0,
                    "final_slope"  : 1.0,
                    "iterations"   : [2,20]
                },
                "change_value_void": {
                    "change_value": false,
                    "iterations"  : [2,100],
                    "initial_value": 0.0,
                    "final_value"  : 0.0
                },
                "change_value_full": {
                    "change_value": false,
                    "iterations"  : [2,100],
                    "initial_value": 1e4,
                    "final_value"  : 1e4
                },
                "domain": ""
            }]
        }""")
        default_physics_parameters_settings.AddMissingParameters(self.GetBasePhysicsParametersSettings())
        return default_physics_parameters_settings
    
    def GetBasePhysicsParametersSettings(self):
        ##settings string in json format
        base_physics_parameters_settings = KratosMultiphysics.Parameters("""
        {
            "resistance": [{}],                                                                
            "conductivity": [{}],
            "decay": [{}],
            "convection_coefficient": [{}],
            "transport_source": [{}]                                                     
        }""")
        return base_physics_parameters_settings
    
    def GetDefaultOptimizationSettings(self):
        ##settings string in json format
        default_optimization_settings = KratosMultiphysics.Parameters("""
        {
            "min_max_iterations": [1,100],
            "change_tolerance"  : 1e-5,
            "optimizer_settings": {
                "optimizer"    : "mma",
                "values_range" : [0.0,1.0],
                "max_outer_it" : 1,
                "kkt_tolerance": 1e-12
            },
            "diffusive_filter_settings": {
                "use_filter" : false,
                "filter_type": "discrete" ,
                "radius"     : 0.01,
                "type_settings": {}
            },
            "projective_filter_settings": {
                "use_filter"              : false,
                "min_max_mean"            : [0.2,0.5],
                "min_max_projection_slope": [1e-10, 2],
                "activation_change"       : 1e-2
            },
            "solution_stabilization_settings": {
                "adjoint_viscosity_adaptation_settings": {
                    "use_adjoint_viscosity_adaptation": false,
                    "adjoint_viscosity_adaptation_factor": 1.0                                                  
                },
                "adjoint_conductivity_adaptation_settings": {
                    "use_adjoint_conductivity_adaptation": false,
                    "adjoint_conductivity_adaptation_factor": 1.0                                                  
                }
            },                                                       
            "remeshing_settings": {
                "remesh"   : false,
                "level_set": 0.15 ,
                "min_max_element_size": [0.01, 0.1]
            },
            "symmetry_settings": {
                "symmetry": false, 
                "model_part_name": ""
            },
            "optimization_domain_settings": 
            {
                "optimization_domain": {
                    "model_part_name": "GENERIC_domain-optimization_domain",
                    "initial_value": 0.0
                },
                "non_optimization_domain": {
                    "model_part_name": "GENERIC_domain-non_optimization_domain",
                    "initial_value": 0.0
                },
                "custom_initial_design_settings": { 
                    "use_custom_initial_design": false,
                    "custom_initial_design": [{
                        "model_part_name": "",
                        "initial_value": "0.0"
                    }]
                }
            },
            "optimization_problem_settings": {
                "functional_weights": {
                    "fluid_functionals": {
                        "time_interval": ["Start","End"],
                        "normalization" : {
                            "type" : "initial",
                            "value": 0.0
                        },
                        "resistance" : {
                            "weight": 1.0,
                            "time_interval": ["Start","End"]
                        },
                        "strain_rate" : {
                            "weight": 1.0,
                            "time_interval": ["Start","End"]
                        },
                        "vorticity" : {
                            "weight": 0.0,
                            "time_interval": ["Start","End"]
                        }
                    },
                    "transport_functionals": {
                        "time_interval": ["Start","End"],
                        "normalization" : {
                            "type" : "initial",
                            "value": 0.0
                        },
                        "outlet_transport_scalar" : {
                            "weight"    : 0.0,
                            "time_interval": ["Start","End"],
                            "outlet_model_part": "Outlet",
                            "target_value": 0.0
                        },
                        "focus_region_transport_scalar" : {
                            "weight"    : 0.0,
                            "time_interval": ["Start","End"],
                            "focus_region_model_part": "FocusRegion",
                            "target_value": 0.0
                        },
                        "conductivity_transfer" : {
                            "weight": 1.0,
                            "time_interval": ["Start","End"]
                        },
                        "convection_transfer" : {
                            "weight": 0.0,
                            "time_interval": ["Start","End"]
                        },
                        "decay_transfer" : {
                            "weight": 0.0,
                            "time_interval": ["Start","End"]
                        },
                        "source_transfer" : {
                            "weight": 0.0,
                            "time_interval": ["Start","End"]
                        },
                        "decay_1st_order_transport" : {
                            "weight": 0.0,
                            "time_interval": ["Start","End"]
                        }
                    },
                    "coupling" : {
                        "fluid"    : 1.0,
                        "transport": 1.0
                    }
                },
                "constraints_settings": {
                    "volume_constraint_settings": {
                        "fluid_or_solid": "fluid",
                        "max_volume_fraction" : 0.4
                    },
                    "use_other_constraints": true,
                    "other_constraints_list": {
                        "WSS_constraint_settings": {
                            "use_WSS_constraint" : false,
                            "min_WSS" : 0.5
                        }
                    }
                }
            },
            "solution_output_settings": {
                "output_path": "solution"
            }
        }""")
        default_optimization_settings.AddMissingParameters(self.GetBaseOptimizationSettings())
        return default_optimization_settings
    
    def GetBaseOptimizationSettings(self):
        ##settings string in json format
        base_physics_parameters_settings = KratosMultiphysics.Parameters("""
        {
            "min_max_iterations": {},
            "change_tolerance"  : {},
            "optimizer_settings": {},
            "diffusive_filter_settings": {},
            "projective_filter_settings": {},
            "solution_stabilization_settings": {},
            "remeshing_settings": {},
            "symmetry_settings": {},
            "optimization_domain_settings": {},
            "optimization_problem_settings": {},
            "solution_output_settings": {}                                                
        }""")
        return base_physics_parameters_settings
    
    def ValidateOptimizationParameters(self):
        """This function validates the settings of the solver
        """
        default_physics_parameters_settings = self.GetDefaultPhysicsParametersSettings()
        self.physics_parameters_settings.ValidateAndAssignDefaults(default_physics_parameters_settings)
        default_optimization_settings = self.GetDefaultOptimizationSettings()
        self.optimization_settings.ValidateAndAssignDefaults(default_optimization_settings)

    def __CheckIfSolveSolutionStepReturnsAValue(self, is_converged):
        """In case the solver does not return the state of convergence
        (same as the SolvingStrategy does) then issue ONCE a deprecation-warning

        """
        if is_converged is None:
            if not hasattr(self, '_map_ret_val_depr_warnings'):
                self._map_ret_val_depr_warnings = []
            solver_class_name = self._GetSolver().__class__.__name__
            # used to only print the deprecation-warning once
            if not solver_class_name in self._map_ret_val_depr_warnings:
                self._map_ret_val_depr_warnings.append(solver_class_name)
                warn_msg  = 'Solver "{}" does not return '.format(solver_class_name)
                warn_msg += 'the state of convergence from "SolveSolutionStep"'
                IssueDeprecationWarning("AnalysisStage", warn_msg)

###########################################################
### METHODS FOR MPI SIMULATION
###########################################################

    def __init__(self,model,parameters):
        self.project_parameters = parameters
        self.topology_optimization_stage = 0
        self.topology_optimization_stage_str = "INIT"
        self._ReadOptimizationParameters()
        super().__init__(model,parameters) 
        self.HandleAdjointSolverAddVariables()
        # self._CreateTopologyOptimizationSolvers() # currently it is a useless method 
        self._SetMinMaxIt()  
        self._SetTopologyOptimizationName()
        self.InitializeDataCommunicator()

    def HandleAdjointSolverAddVariables(self):
        self._GetAdjointSolver().AddVariables()
        KratosMultiphysics.MergeVariableListsUtility().Merge(self._GetPhysicsSolver().main_model_part, self._GetAdjointSolver().main_model_part)

    def InitializeDataCommunicator(self):
        self.data_communicator = DataCommunicator.GetDefault()

    def _CreateNodesIdsDictionary(self):
        """
        Creates dictionaries mapping between global and local node IDs for the current MPI partition.
        The dictionary `nodes_ids_global_to_local_partition_dictionary` maps global node IDs to local indices,
        while `nodes_ids_local_partition_to_global_dictionary` provides the reverse mapping.
        These mappings are used to ensure consistent node indexing across MPI partitions.
        If the simulation is not running in MPI, the global and local IDs coincide,
        and the dictionaries define an identity mapping.
        """
        self.MpiPrint("--|" + self.topology_optimization_stage_str + "| ---> Creates Nodes Ids Dictionary")
        self.nodes_ids_global_to_local_partition_dictionary = {}
        self.nodes_ids_local_partition_to_global_dictionary = {}
        count = 0
        for node in self._GetLocalMeshNodes():
            self.nodes_ids_global_to_local_partition_dictionary[node.Id] = count
            self.nodes_ids_local_partition_to_global_dictionary[count]   = node.Id
            count += 1
        if (count != len(self.nodes_ids_global_to_local_partition_dictionary.keys())):
            raise RuntimeError("Wrong reordering of nodes ids. The counted number of nodes is different from len(self._GetLocalMeshNodes()).")
        self._PassNodesIdsGlobalToLocalDictionaryToSolvers()
        self.n_nodes = count
        self._EvaluateTotalNumberOfNodes()

    def _PassNodesIdsGlobalToLocalDictionaryToSolvers(self):
        self._GetPhysicsSolver().SetNodesIdsGlobalToLocalDictionary(self.nodes_ids_global_to_local_partition_dictionary)
        self._GetAdjointSolver().SetNodesIdsGlobalToLocalDictionary(self.nodes_ids_global_to_local_partition_dictionary)

    def _EvaluateTotalNumberOfNodes(self):
        local_n_nodes = self.n_nodes
        total_n_nodes = self.data_communicator.SumAll(local_n_nodes)
        self.total_n_nodes = total_n_nodes

    def _CreateElementsIdsDictionary(self):
        self.MpiPrint("--|" + self.topology_optimization_stage_str + "| ---> Creates Elements Ids Dictionary")
        self.elements_ids_dictionary = {}
        count = 0
        mp = self._GetComputingModelPart()
        for element in mp.GetCommunicator().LocalMesh().Elements:
            self.elements_ids_dictionary[element.Id] = count
            count += 1
        if (count != len(self.elements_ids_dictionary.keys())):
            raise RuntimeError("Wrong reordering of nodes ids. The counted number of nodes is different from len(mp.GetCommunicator().LocalMesh().Elements).")
        self.n_elements = count
        self._EvaluateTotalNumberOfElements()

    def _EvaluateTotalNumberOfElements(self):
        local_n_elements = self.n_elements
        total_n_elements = self.data_communicator.SumAll(local_n_elements)
        self.total_n_elements = total_n_elements

    def _EvaluateFunctional(self, print_functional=False):
        """
        This method is used to evaluate the functional value
        # Functionals Database
        # 0: resistance  : int_{\Omega}{alpha*||u||^2}
        # 1: strain-rate : int_{\Omega}{2*mu*||S||^2} , with S = 1/2*(grad(u)+grad(u)^T) strain-rate tensor
        # 2: vorticity   : int_{\Omega}{2*mu*||R||^2} = int_{\Omega}{mu*||curl(u)||^2} , curl(u) = vorticity vector, R = 1/2*(grad(u)-grad(u)^T) rotation-rate tensor
        # 3: outlet_transport_scalar : int_{\Gamma_{out}}{c}
        # 4: region_transport_scalar: int_{\Omega}{c^2}
        # 5: transport_scalar_diffusion: int_{\Omega}{D\\||grad(u)||^2}
	    # 6: transport_scalar_convection: int_{\Omega}{beta*T*dot(u,grad(T))}
	    # 7: transport_scalar_decay: int_{\Omega}{kT^2}
	    # 8: transport_scalar_source: int_{\Omega}{-Q*T}
        """
        self._SetTopologyOptimizationStage(3)
        self.MpiPrint("--|" + self.topology_optimization_stage_str + "| EVALUATE FUNCTIONAL VALUE")
        self.EvaluateFunctionals(print_functional)
        self.EvaluateTotalFunctional()

    def _EvaluateVolumeConstraintAndDerivative(self):
        self.EvaluateDesignParameterIntegralInOptimizationDomain()
        if (self.is_fluid_volume_constraint):
            self.volume_fraction = 1.0 - self.design_parameter_integral/self.optimization_domain_size
            volume_constraint_derivatives_wrt_design_base = -1.0 * self.nodal_optimization_domain_sizes / self.optimization_domain_size
        else:
            self.volume_fraction = self.design_parameter_integral/self.optimization_domain_size
            volume_constraint_derivatives_wrt_design_base = self.nodal_optimization_domain_sizes / self.optimization_domain_size
        self.volume_constraint = self.volume_fraction - self.max_volume_fraction
        volume_constraint_derivatives_wrt_design_projected = volume_constraint_derivatives_wrt_design_base * self.design_parameter_projected_derivatives
        self.constraints[self.volume_constraint_id] = self.volume_constraint
        self.constraints_derivatives_wrt_design[self.volume_constraint_id,:] = self._ApplyDiffusiveFilterDerivative(volume_constraint_derivatives_wrt_design_projected)

    def EvaluateDesignParameterIntegralInOptimizationDomain(self):
        if self.IsMpiParallelism():
            local_design_parameter_integral = np.dot(self.design_parameter, self.nodal_optimization_domain_sizes)
            total_design_parameter_integral = self.data_communicator.SumAll(local_design_parameter_integral)
            self.design_parameter_integral  = total_design_parameter_integral
        else:
            self.design_parameter_integral  = np.dot(self.design_parameter, self.nodal_optimization_domain_sizes)

    def _ComputeDomainSize(self):
        self.MpiPrint("--|" + self.topology_optimization_stage_str + "| ---> Compute Domain Size")
        self._ComputeOptimizationDomainSize() # must be done before since it uses the compute nodal area process that after has to be overwritten by the one evaluated in the whole domain model part
        self._ComputeElementalDomainSize()
        self._ComputeNodalDomainSize()

    def _ComputeOptimizationDomainSize(self):
        self._InitializeOptimizationDomainSettings()
        self._ComputeOptimizationDomainNodesMask()
        self._ComputeElementalOptimizationDomainSize()
        self._ComputeNodalOptimizationDomainSize() 

    def _ComputeElementalDomainSize(self):
        self.elemental_domain_size = np.zeros(self.n_elements)
        self.total_domain_size = 0.0
        mp = self._GetComputingModelPart()
        for elem in mp.Elements:
            temp_domain_size = elem.GetGeometry().DomainSize()
            self.elemental_domain_size[self.elements_ids_dictionary[elem.Id]] = temp_domain_size
        local_domain_size = np.sum(self.elemental_domain_size)
        # synchronize the value across all the ranks of the data_communicator
        # self.MpiBarrier()
        total_area = mp.GetCommunicator().GetDataCommunicator().SumAll(local_domain_size)
        self.total_domain_size = total_area

    def _ComputeElementalOptimizationDomainSize(self):
        self.optimization_domain_size = 0.0
        local_domain_size = 0.0
        mp = self._GetOptimizationDomain()
        for elem in mp.Elements:
            local_domain_size += elem.GetGeometry().DomainSize()
        # synchronize the value across all the ranks of the data_communicator
        # self.MpiBarrier()
        total_domain_size = mp.GetCommunicator().GetDataCommunicator().SumAll(local_domain_size)
        self.optimization_domain_size = total_domain_size

    def _ComputeNodalDomainSize(self):
        self.nodal_domain_sizes = np.zeros(self.n_nodes)
        mp = self._GetComputingModelPart()
        nodal_area_process = KratosMultiphysics.CalculateNodalAreaProcess(mp, self.dim)
        nodal_area_process.Execute()
        self._UpdateNodalDomainSizeArrayFromNodalAreaVariable()
        self._CorrectNodalDomainSizeWithSymmetry()

    def _ComputeNodalOptimizationDomainSize(self):
        self.nodal_optimization_domain_sizes = np.zeros(self.n_nodes)
        mp = self._GetOptimizationDomain()
        if self.CurrentDomainHasOptimizationNodes():
            nodal_area_process = KratosMultiphysics.CalculateNodalAreaProcess(mp, self.dim)
            nodal_area_process.Execute()
            self._UpdateNodalOptimizationDomainSizeArrayFromNodalAreaVariable()
            self._CorrectNodalOptimizationDomainSizeWithSymmetry()

    def _UpdateNodalDomainSizeArrayFromNodalAreaVariable(self):
        for node in self._GetLocalMeshNodes():
            self.nodal_domain_sizes[self.nodes_ids_global_to_local_partition_dictionary[node.Id]] = node.GetSolutionStepValue(KratosMultiphysics.NODAL_AREA)

    def _UpdateNodalOptimizationDomainSizeArrayFromNodalAreaVariable(self):
        for node in self._GetLocalMeshNodes(self._GetOptimizationDomain()):
            self.nodal_optimization_domain_sizes[self.nodes_ids_global_to_local_partition_dictionary[node.Id]] = node.GetSolutionStepValue(KratosMultiphysics.NODAL_AREA)

    def _ComputeOptimizationDomainNodesMask(self):
        """
        Builds the mapping (mask) between the total set of domain nodes and the subset of nodes
        belonging to the optimization domain.
        The method creates two arrays:
        - `global_to_opt_mask_nodes`: an array of size `n_nodes` that maps each node index
        to its position in the optimization-domain mask, or to -1 if the node belongs to the
        non-optimization domain.
        - `optimization_domain_nodes_mask`: an array containing the indices of the nodes
        that belong to the optimization domain only.
        This mask is used to transfer quantities between the full computational domain and
        the optimization domain.
        """
        self.MpiPrint("--|" + self.topology_optimization_stage_str + "| ---> Compute Optimization Domain Nodes Mask")
        non_opt_mp= self._GetNonOptimizationDomain()
        # opt_mp= self._GetOptimizationDomain()
        mp = self._GetComputingModelPart()
        #define the boolean search for the nodes belonging only to the optimization domain model part
        self.global_to_opt_mask_nodes = np.ones(self.n_nodes, dtype=int)
        # -1: non opt
        # else: position in the optimization mask
        if (not (non_opt_mp is None)):
            for node in self._GetLocalMeshNodes(non_opt_mp):
                self.global_to_opt_mask_nodes[self.nodes_ids_global_to_local_partition_dictionary[node.Id]] = -1
            n_non_opt_nodes = len(self._GetLocalMeshNodes(non_opt_mp))
        else:
            n_non_opt_nodes = 0
        self.n_opt_design_parameters = self.n_nodes - n_non_opt_nodes
        self.optimization_domain_nodes_mask = np.zeros(self.n_opt_design_parameters, dtype=int)
        count_opt_nodes = 0
        for node in self._GetLocalMeshNodes(mp):
            full_domain_to_local_rank_node_id = self.nodes_ids_global_to_local_partition_dictionary[node.Id]
            if (self.global_to_opt_mask_nodes[full_domain_to_local_rank_node_id] != -1): #if the node is only in the optimization domain, add it to the mask
                self.optimization_domain_nodes_mask[count_opt_nodes] = full_domain_to_local_rank_node_id
                self.global_to_opt_mask_nodes[full_domain_to_local_rank_node_id] = count_opt_nodes
                count_opt_nodes +=1
        if (count_opt_nodes != self.n_opt_design_parameters):
            self.MpiPrint("!!! WARNING: wrong initialization of the Optimization Domain Nodes Mask")

    def _UpdateFunctionalDerivativesVariable(self):
        for node in self._GetLocalMeshNodes():
            node.SetValue(KratosMultiphysics.FUNCTIONAL_DERIVATIVE, self.functional_derivatives_wrt_design[self.nodes_ids_global_to_local_partition_dictionary[node.Id]][0])

    def _UpdateConstraintsDerivativesVariable(self):
        self._UpdateVolumeConstraintDerivativesVariable()
    
    def _UpdateVolumeConstraintDerivativesVariable(self):
        # for node in self._GetLocalMeshNodes():
        #     node.SetValue(KratosMultiphysics.FUNCTIONAL_DERIVATIVE, self.functional_derivatives_wrt_design[self.nodes_ids_global_to_local_partition_dictionary[node.Id]][0])
        pass
    
    def MpiSumLocalValues(self, local_value):
        """
        Performs a global sum reduction of the provided local value across all MPI ranks.
        The method uses `data_communicator.SumAll()` to compute the sum of `local_value`
        over all partitions and returns the total value.
        It is used to pass to every rank the total value of the functional.
        """
        total_value = self.data_communicator.SumAll(local_value)
        return total_value

    def UpdatePhysicsParametersVariablesAndSynchronize(self):
        self.UpdatePhysicsParametersVariables()
        self._SynchronizePhysicsParametersVariables()
        
    def _SynchronizePhysicsParametersVariables(self):
        self._GetMainModelPart().GetCommunicator().SynchronizeNonHistoricalVariable(KratosCFD.RESISTANCE)

    def _SolveMMA(self, design_parameter, n_opt_variables, n_opt_constraints, min_value, max_value, max_outer_it, kkt_tolerance):
        """
        Solves the MMA (Method of Moving Asymptotes) subproblem in parallel.
        
        Algorithm:
        - Initializes MMA parameters and variable bounds (xmin/xmax, low/upp, move, a0, a, c, d).
        - Iterates the MMA outer loop until the KKT norm is below `kkt_tolerance` or the maximum number of iterations
        (`max_outer_it`) is reached:
            * Solves the MMA subproblem using `MMA.mmasub(...)` to update `xval`.
            * Re-evaluates the objective and constraints and their gradients.
            * Checks the KKT residuals with `MMA.kktcheck(...)`.
        - Extracts the portion of `xval` corresponding to the current rank and maps it back to the local design parameter
        vector using `_InsertDesignParameterFromOptimizationDomain(...)`, then updates the model and physical parameters
        through `_UpdateDesignParameterAndPhysicsParameters(...)`.

        MPI logic:
        - Gathers per-rank sizes with `AllGatherInts([n_opt_variables])` to obtain the global number of design variables
        and the rank offsets.
        - Gathers the per-rank design vectors with `AllGathervDoubles(design_parameter)` and concatenates them into the
        global design vector `xval` of size N = sum(n_opt_variables_in_rank).
        - Each rank operates on its own slice `xval[offset[rank]:offset[rank+1]]` when evaluating the objective and
        constraint functions via `_UpdateOptimizationProblem`, while the corresponding gradients are assembled globally
        with `CreateFunctionalAndConstraintsDerivativesCompleteArrays(...)`.
        As a result, all ranks hold identical global arrays for `f0val`, `df0dx`, `fval`, and `dfdx`.

        - In serial execution, the gather operations are trivial: `N == n_opt_variables` and the offsets are [0, N].
        - `min_value` and `max_value` define the global bounds applied to all design variables.
        """
        self.MpiPrint("--|" + self.topology_optimization_stage_str + "| SOLVE MMA")
        # MMA PARAMETERS INITIALIZATION
        rank = self.data_communicator.Rank()
        n_opt_variables_in_rank = self.data_communicator.AllGatherInts([n_opt_variables])
        n_ranks = len(n_opt_variables_in_rank)
        total_opt_variables = sum(n_opt_variables_in_rank)
        design_parameter_per_rank = self.data_communicator.AllGathervDoubles(design_parameter)
        get_current_rank_nodes_subset = [0]
        for irank in range(len(n_opt_variables_in_rank)):
            get_current_rank_nodes_subset.append(get_current_rank_nodes_subset[irank] + n_opt_variables_in_rank[irank])
        # self.MpiBarrier()
        xval = np.concatenate(design_parameter_per_rank)
        # self.MpiBarrier()
        n = total_opt_variables
        m = n_opt_constraints
        eeen = np.ones((n, 1))
        eeem = np.ones((m, 1))
        zeron = np.zeros((n, 1))
        zerom = np.zeros((m, 1))
        xval  = xval.reshape(-1, 1)
        xold1 = xval.copy()
        xold2 = xval.copy()
        xmin = min_value * eeen
        xmax = max_value * eeen
        low = xmin.copy()
        upp = xmax.copy()
        move = 0.4
        c = 1000 * eeem
        d = eeem.copy()
        a0 = 1
        a = zerom.copy()
        innerit = 0
        outeriter = 0
        maxoutit = max_outer_it
        kkttol = kkt_tolerance
        # Calculate function values and gradients of the objective and constraints functions
        if (outeriter == 0):
            # self.MpiBarrier()
            curr_rank_design = xval[get_current_rank_nodes_subset[rank]:get_current_rank_nodes_subset[rank+1]]
            f0val, temp_df0dx, fval, temp_dfdx = self._UpdateOptimizationProblem(curr_rank_design.flatten())
            df0dx, dfdx = self.CreateFunctionalAndConstraintsDerivativesCompleteArrays(temp_df0dx, temp_dfdx, n_ranks, n_opt_variables_in_rank, n_opt_constraints)
            # self.MpiBarrier()
        # The iterations start
        kktnorm = kkttol + 10
        outit = 0
        while ((kktnorm > kkttol) and (outit < maxoutit)):
            outit += 1
            outeriter += 1
            # The MMA subproblem is solved at the point xval:
            xmma, ymma, zmma, lam, xsi, eta, mu, zet, s, low, upp = MMA.mmasub(
                m, n, outeriter, xval, xmin, xmax, xold1, xold2, f0val, df0dx, fval, dfdx, low, upp, a0, a, c, d, move)
            # Store previous results:
            xold2 = xold1.copy()
            xold1 = xval.copy()
            xval = xmma.copy()
            # Re-calculate function values and gradients of the objective and constraints functions
            # self.MpiBarrier()
            curr_rank_design = xval[get_current_rank_nodes_subset[rank]:get_current_rank_nodes_subset[rank+1]]
            f0val, temp_df0dx, fval, temp_dfdx = self._UpdateOptimizationProblem(curr_rank_design.flatten())
            df0dx, dfdx = self.CreateFunctionalAndConstraintsDerivativesCompleteArrays(temp_df0dx, temp_dfdx, n_ranks, n_opt_variables_in_rank, n_opt_constraints)
            # self.MpiBarrier()
            # The residual vector of the KKT conditions is calculated
            residu, kktnorm, residumax = MMA.kktcheck(
                m, n, xmma, ymma, zmma, lam, xsi, eta, mu, zet, s, xmin, xmax, df0dx, fval, dfdx, a0, a, c, d)
        curr_rank_design = xval[get_current_rank_nodes_subset[rank]:get_current_rank_nodes_subset[rank+1]]
        new_design_parameter = self._InsertDesignParameterFromOptimizationDomain(curr_rank_design.flatten())
        self._UpdateDesignParameterAndPhysicsParameters(new_design_parameter)

    def CreateFunctionalDerivativesCompleteArrays(self, local_value, n_ranks, n_opt_variables_in_rank):
        """
        Gathers and concatenates the objective function derivatives (df0/dx) from all MPI ranks
        into a single global array.

        Each rank provides its local gradient `local_value`, which is gathered using
        `AllGathervDoubles`. The gathered arrays are then reshaped according to the number
        of optimization variables per rank (`n_opt_variables_in_rank`) and concatenated
        along the first axis to form the complete global derivative vector.

        In serial execution, this simply returns `local_value` as a column vector.
        """
        local_value = local_value.flatten().tolist()
        gather_local_value = self.data_communicator.AllGathervDoubles(local_value)
        for irank in range(n_ranks):
            gather_local_value[irank] = np.asarray(gather_local_value[irank]).reshape(n_opt_variables_in_rank[irank],1)
        result = np.concatenate(gather_local_value, axis=0)
        return result
    
    def CreateConstraintsDerivativesCompleteArrays(self, local_value, n_ranks, n_opt_variables_in_rank, n_opt_constraints):
        """
        Gathers and concatenates the constraint function derivatives (df/dx) from all MPI ranks
        into a single global array.

        Each rank provides its local derivative matrix `local_value`, which is gathered using
        `AllGathervDoubles`. The gathered matrices are reshaped to size
        (`n_opt_constraints`, `n_opt_variables_in_rank[irank]`) and concatenated along the
        second axis to form the complete global constraint Jacobian.

        In serial execution, this simply returns `local_value` reshaped to
        (`n_opt_constraints`, `n_opt_variables_in_rank[0]`).
        """
        local_value = local_value.flatten().tolist()
        gather_local_value = self.data_communicator.AllGathervDoubles(local_value)
        for irank in range(n_ranks):
            gather_local_value[irank] = np.asarray(gather_local_value[irank]).reshape(n_opt_constraints, n_opt_variables_in_rank[irank])
        result = np.concatenate(gather_local_value, axis=1)
        return result
    
    def CreateFunctionalAndConstraintsDerivativesCompleteArrays(self, temp_df0dx, temp_dfdx, n_ranks, n_opt_variables_in_rank, n_opt_constraints):
        """
        Builds the complete global derivative arrays for both the objective function and
        the constraints by combining the contributions from all MPI ranks.
        - Calls `CreateFunctionalDerivativesCompleteArrays` to assemble the global objective
        gradient `df0dx`.
        - Calls `CreateConstraintsDerivativesCompleteArrays` to assemble the global constraint
        Jacobian `dfdx`.
        Returns:
            df0dx (np.ndarray): Global gradient of the objective function (shape: [N, 1]).
            dfdx (np.ndarray): Global derivatives of the constraints (shape: [M, N]).
        In serial execution, the returned arrays coincide with the local derivatives.
        """
        df0dx = self.CreateFunctionalDerivativesCompleteArrays(temp_df0dx, n_ranks, n_opt_variables_in_rank)
        dfdx = self.CreateConstraintsDerivativesCompleteArrays(temp_dfdx, n_ranks, n_opt_variables_in_rank, n_opt_constraints)
        return df0dx, dfdx
    
    def _GetViscosity(self):
        """
        Returns the dynamic viscosity of the fluid from the main model part.
        The method accesses the first element in the local mesh and retrieves the value
        of `DYNAMIC_VISCOSITY` from its properties, assuming a uniform viscosity across
        all elements.
        """
        for elem in self._GetMainModelPart().GetCommunicator().LocalMesh().Elements:
            mu = elem.Properties.GetValue(KratosMultiphysics.DYNAMIC_VISCOSITY)
            break
        return mu
    
    def _GetDensity(self):
        """
        Returns the density of the fluid from the main model part.

        The method accesses the first element in the local mesh and retrieves the value
        of `DENSITY` from its properties, assuming a uniform density across all elements.
        """
        for elem in self._GetMainModelPart().GetCommunicator().LocalMesh().Elements:
            rho = elem.Properties.GetValue(KratosMultiphysics.DENSITY)
            break
        return rho
    
    def _GetLocalMeshNodes(self, mp = None):
        if mp is None:
            return self._GetPhysicsMainModelPart().GetCommunicator().LocalMesh().Nodes
        else:
            return mp.GetCommunicator().LocalMesh().Nodes
        
    def _GetLocalMeshElements(self, mp = None):
        if mp is None:
            return self._GetMainModelPart().GetCommunicator().LocalMesh().Elements
        else:
            return mp.GetCommunicator().LocalMesh().Elements
        
    def _GetPhysicsSolverDistributedModelPartImporter(self):
        return self.physics_solver.distributed_model_part_importer

    def _CorrectPvtuFilesInOptimizationVtuOutput(self):
        """
        Corrects the relative paths in the .pvtu files generated during VTU output.
        This method should be executed only on rank 0. It checks the `vtu_output`
        settings in `project_parameters` and identifies the `.pvtu` file corresponding
        to the current optimization iteration (`self.opt_it`).
        If the file exists, it reads its content and removes redundant folder path
        prefixes (e.g., "output_path/") from the internal VTU file references, ensuring
        that the `.pvtu` file correctly links to the corresponding `.vtu` pieces.
        This fix is necessary because the default Kratos VTU output may include
        incorrect relative paths when writing parallel `.pvtu` files.
        """
        if (self.MpiRunOnlyRank(0)):
            output_processes_list_names = ["optimization_output_processes"]
            for output_process_list_name in output_processes_list_names:
                if (self.project_parameters.Has(output_process_list_name)):
                    if (self.project_parameters[output_process_list_name].Has("vtu_output")):
                        vtu_output_settings = self.project_parameters[output_process_list_name]["vtu_output"][0]
                        if (vtu_output_settings.Has("Parameters")):
                            vtu_output_parameters = vtu_output_settings["Parameters"]
                            mp_name     = vtu_output_parameters["model_part_name"].GetString()
                            folder_name = vtu_output_parameters["output_path"].GetString()
                            curr_it_str = str(self.opt_it)
                            file_name = f"{mp_name}_{curr_it_str}.pvtu"
                            pvtu_path = Path(folder_name) / file_name
                            if pvtu_path.exists():
                                text = pvtu_path.read_text(encoding='utf-8')
                                fixed_text = text.replace(folder_name+"/", "")
                                pvtu_path.write_text(fixed_text, encoding='utf-8')

    def _CorrectPvtuFilesInTimeVtuOutput(self):
        if (self.MpiRunOnlyRank(0)):
            output_processes_list_names = ["time_output_processes"]
            for output_process_list_name in output_processes_list_names:
                if (self.project_parameters.Has(output_process_list_name)):
                    if (self.project_parameters[output_process_list_name].Has("vtu_output")):
                        vtu_output_settings = self.project_parameters[output_process_list_name]["vtu_output"][0]
                        if (vtu_output_settings.Has("Parameters")):
                            vtu_output_parameters = vtu_output_settings["Parameters"]
                            folder_name = vtu_output_parameters["output_path"].GetString()
                            for istep in range(self.n_time_steps):
                                time_step_counter = str(istep+1).zfill(len(str(self.n_time_steps))) 
                                file_name = "time_step_" + time_step_counter + ".pvtu"
                                pvtu_path = Path(folder_name) / file_name
                                if pvtu_path.exists():
                                    text = pvtu_path.read_text(encoding='utf-8')
                                    fixed_text = text.replace(folder_name+"/", "")
                                    pvtu_path.write_text(fixed_text, encoding='utf-8')

###########################################################
### METHODS FOR MPI UTILITIES
###########################################################

    def MpiCheck(self, text="before solving step", rank=-1):
        if (rank == -1): # print for all ranks
            print("--|" + str(self.data_communicator.Rank()) + "| Checkpoint reached:", text)
        elif (self.data_communicator.Rank() == rank): # print only for a specific rank
            print("--|" + str(rank) + "| Checkpoint reached:", text)

    def MpiBarrier(self):
        self.data_communicator.Barrier()

    def MpiCheckAndBarrier(self, text="before solving step", rank=-1):
        # self.MpiBarrier()
        self.MpiCheck(text=text, rank=rank)
        # self.MpiBarrier()
    
    def MpiRunOnlyRank(self, rank=0):
        """
        Returns: True if the simulation is not distributed or if it is running on a specified data_communicator rank
        """
        if (self.IsOpenMPParallelism()):
            return True
        elif (self.data_communicator.Rank() == rank):
            return True
        else:
            return False
        
    def IsMpiParallelism(self):
        return _CheckIsDistributed()
    
    def IsOpenMPParallelism(self):
        return not self.IsMpiParallelism()
    
    def MpiNotRunOnlyRank(self, rank=0):
        return not self.MpiRunOnlyRank(rank)
        
    def MpiPrint(self, text_to_print="", rank=0, set_barrier=False, min_echo=1):
        if self.echo_level >= min_echo:
            if (self.IsOpenMPParallelism()):
                print(text_to_print)
            else:
                if (set_barrier):
                    self.MpiBarrier()
                if (self.MpiRunOnlyRank(rank)):
                    print(text_to_print)
                if (set_barrier):
                    self.MpiBarrier()   
        else:
            pass

    def CurrentDomainHasOptimizationNodes(self):
        opt_mp = self._GetOptimizationDomain()
        if (len(opt_mp.Nodes) != 0):
            return True
        else:
            return False
        
    def GetPhysicsBufferIdFromTimeStepId(self, time_step_id):
        # assume time step id starts from 0 for the first time step and ends with self.n_time_steps-1
        return self.n_time_steps-(time_step_id+1)
    
    def GetAdjointBufferIdFromTimeStepId(self, time_step_id):
        # assume time step id starts from 0 for the first time step and ends with self.n_time_steps-1
        return time_step_id
    
    def _SolveFluidPhysics(self):
        return True
    
    def _SolveTransportPhysics(self):
        return False
