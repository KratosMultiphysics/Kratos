## IMPORT
# Import Libraries
from sys import argv
import numpy as np #import the numpy library
import scipy as sp #import the scipy library
from scipy.spatial import KDTree
from scipy.sparse import dok_matrix, lil_matrix
import mmapy as MMA # import the MMA subroutines python library: https://github.com/arjendeetman/GCMMA-MMA-Python
import time

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
from KratosMultiphysics.FluidDynamicsApplication import fluid_topology_optimization_solver
from KratosMultiphysics.ConvectionDiffusionApplication.apply_topology_optimization_pde_filter_process_mpi import ApplyTopologyOptimizationPdeFilterProcessMpi
# Import Kratos Processes
from KratosMultiphysics import ComputeNodalGradientProcess

class FluidTopologyOptimizationAnalysisMpi(FluidDynamicsAnalysis):
    def _ReadOptimizationParameters(self):
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
        self.physics_solver = fluid_topology_optimization_solver.CreateSolver(self.model, self.project_parameters)
        self.adjoint_solver = fluid_topology_optimization_solver.CreateSolver(self.model, self.project_parameters, isAdjointSolver=True)

    def _CreateSolver(self, isAdjointSolver = False):
        """
        This method creates a solver
        isAdjointSolver == False --> physics_solver
        isAdjointSolver == True  --> adjoint_solver
        """
        return fluid_topology_optimization_solver.CreateSolver(self.model, self.project_parameters, isAdjointSolver)
    
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
        
    def _SetMinMaxIt(self):
        """
        This method sets the minimum & maximum iteration for the topology optimization process
        """
        iterations_settings = self.optimization_settings["min_max_iterations"]
        self.min_it = iterations_settings[0].GetInt()
        self.max_it = iterations_settings[1].GetInt()
        if (self.min_it > self.max_it):
            self.MpiPrint("\n!!!WARNING: wrong initialization of the min & max number of iterations\n")

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

    def _SetFunctionalWeights(self):
        # set future transport functinals to zero
        self._InitializeFunctionalWeights()
        self._GetComputingModelPart().ProcessInfo.SetValue(KratosMultiphysics.FUNCTIONAL_WEIGHTS, self.functional_weights)
        self._PrintFunctionalWeights()

    def _PrintFunctionalWeights(self):
        self.MpiPrint("--|" + self.topology_optimization_stage_str + "| FUNCTIONAL WEIGHTS: " + str(self.functional_weights))
        self.MpiPrint(self.optimization_settings["optimization_problem_settings"]["functional_weights"].PrettyPrintJsonString())
        
    def _InitializeFunctionalWeights(self):
        fluid_weights = self._ImportFluidFunctionalWeights()
        weights = fluid_weights + [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
        self.n_functionals = len(weights)
        self.initial_functionals_values = np.zeros(self.n_functionals)
        self.functionals = np.zeros(self.n_functionals)
        self.functional_weights = self._NormalizeFunctionalWeights(np.asarray(weights))

    def _ImportFluidFunctionalWeights(self):
        fluid_weights = [0.0, 0.0, 0.0]
        functional_weights_parameters = self.optimization_settings["optimization_problem_settings"]["functional_weights"]["fluid_functionals"]
        fluid_weights[0] = functional_weights_parameters["resistance"]["weight"].GetDouble()
        fluid_weights[1] = functional_weights_parameters["strain_rate"]["weight"].GetDouble()
        fluid_weights[2] = functional_weights_parameters["vorticity"]["weight"].GetDouble()
        return fluid_weights
    

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
        self._InitializeTopologyOptimizationProblem()
        self.SolveTopologyOptimization()

    def _InitializeTopologyOptimizationProblem(self):
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
        # self._ResetFunctionalOutput()

    def SolveTopologyOptimization(self):
        self.design_parameter_change = self.design_parameter_change_toll + 10
        end_solution = False
        self.first_iteration = True
        while (not end_solution):
            self.opt_it = self.opt_it+1
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
            self._PrintSolution()
            self.first_iteration = False
            
    def _IsTopologyOptimizationSolutionEnd(self):
        design_parameter_converged = self._EvaluateDesignParameterChange()
        volume_constraint_valid = not (self.volume_constraint > 0.0)
        self.converged = design_parameter_converged and volume_constraint_valid
        return not (((self.opt_it < self.max_it) and (not self.converged)) or (self.opt_it < self.min_it))

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
            self.OutputSolutionStep()
    
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
        self._ComputeDesignParameterFilterUtilities()
        self.opt_it = 0
        self._SetDesignParameterChangeTolerance()
        self.n_optimization_constraints = 0  
        self._InitializeOptimizerSettings()
        self._InitializeConstraints()
        self._InitializeRemeshing()
    
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
        mp = self._GetComputingModelPart()
        counter = 0
        for node in mp.GetCommunicator().LocalMesh().Nodes:
            design = self.design_parameter[counter]
            counter += 1
            node.SetSolutionStepValue(KratosMultiphysics.DESIGN_PARAMETER, design)
            distance = design-self.remeshing_levelset
            node.SetSolutionStepValue(KratosMultiphysics.DISTANCE, distance)

    def _SetDesignParameterCustomInitialDesign(self):
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
                self.symmetry_nodes_mask = np.zeros(len(sym_mp.GetCommunicator().LocalMesh().Nodes), dtype=int)
                count = 0
                for node in sym_mp.GetCommunicator().LocalMesh().Nodes:
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

    def _InitializeOptimizationDomainSettings(self):
        optimization_domain_settings = self.optimization_settings["optimization_domain_settings"]
        self.optimization_domain_name = optimization_domain_settings["optimization_domain"]["model_part_name"].GetString()
        self.optimization_domain_initial_value = max(0.0, min(1.0, optimization_domain_settings["optimization_domain"]["initial_value"].GetDouble()))
        self.non_optimization_domain_name = optimization_domain_settings["non_optimization_domain"]["model_part_name"].GetString()
        self.non_optimization_domain_initial_value = max(0.0, min(1.0, optimization_domain_settings["non_optimization_domain"]["initial_value"].GetDouble()))

    def _CorrectNodalOptimizationDomainSizeWithSymmetry(self):
        if (self.symmetry_enabled):
            self.nodal_optimization_domain_sizes[self.symmetry_nodes_mask] *= 2.0

    def _GetModelPartNodesSubset(self, model_part, node_ids):
        return [model_part.GetNode(node_id) for node_id in node_ids]
    
    def _GetModelPartNodesIds(self, model_part):
        return [node.Id for node in model_part.Nodes]

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
        self.shape_functions_derivatives = np.zeros((self.n_elements, self.nodes_in_element, self.dim))
        self.element_nodes_ids = np.zeros((self.n_elements, self.nodes_in_element), dtype=int)
        mp = self._GetComputingModelPart()
        for el in mp.Elements:
            el_geometry = el.GetGeometry()
            global_gradients_at_nodes = self._GetShapeFunctionsDerivatives(el)
            count_node = 0
            for node in el_geometry:
                self.shape_functions_derivatives[self.elements_ids_dictionary[el.Id], count_node, :] = global_gradients_at_nodes[count_node,:]
                self.element_nodes_ids[self.elements_ids_dictionary[el.Id], count_node] = self.nodes_ids_global_to_local_partition_dictionary[node.Id]
                count_node = count_node+1

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
        mp = self._GetComputingModelPart()
        for node in mp.GetCommunicator().LocalMesh().Nodes:
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
        domain = physics_parameters["domain"].GetString()
        mp = self._GetSubModelPart(self._GetMainModelPart(), domain)
        if (mp is not None):
                nodes_ids = self.nodes_ids_global_to_local_partition_dictionary[self._GetModelPartNodesIds(mp)]
        else:
            nodes_ids = np.arange(self.n_nodes)
        interpolation_method = (physics_parameters["interpolation_method"].GetString()).lower()
        if (interpolation_method == "hyperbolic"):
            computed_parameter, computed_parameter_derivative_wrt_design_base = self._ComputeHyperbolicPhysicsParameter(physics_parameters, design_parameter)
        elif (interpolation_method == "polynomial"):
            computed_parameter, computed_parameter_derivative_wrt_design_base = self._ComputePolynomialPhysicsParameter(physics_parameters, design_parameter)
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
        # Adjust void value
        if (physics_parameters["change_value_void"]["change_value"].GetBool()):
            start_it      = physics_parameters["change_value_void"]["iterations"][0].GetInt()
            end_it      = physics_parameters["change_value_void"]["iterations"][1].GetInt()
            initial_value = physics_parameters["change_value_void"]["initial_value"].GetDouble()
            final_value = physics_parameters["change_value_void"]["final_value"].GetDouble()
            if (self.opt_it < self.min_it):
                current_value = initial_value
            elif (self.opt_it < end_it):
                current_value = initial_value + (final_value-initial_value)*(self.opt_it-start_it)/(end_it-start_it)
            else:
                current_value = final_value
            physics_parameters["value_void"].SetDouble(current_value)

        # Adjust full value
        if (physics_parameters["change_value_full"]["change_value"].GetBool()):
            start_it      = physics_parameters["change_value_full"]["iterations"][0].GetInt()
            end_it      = physics_parameters["change_value_full"]["iterations"][1].GetInt()
            initial_value = physics_parameters["change_value_full"]["initial_value"].GetDouble()
            final_value = physics_parameters["change_value_full"]["final_value"].GetDouble()
            if (self.opt_it < start_it):
                current_value = initial_value
            elif (self.opt_it < end_it):
                current_value = initial_value + (final_value-initial_value)*(self.opt_it-start_it)/(end_it-start_it)
            else:
                current_value = final_value
            physics_parameters["value_full"].SetDouble(current_value)



    
    def _UpdateResistanceDesignDerivative(self):
        mask = self._GetOptimizationDomainNodesMask()
        resistance_derivative_wrt_design_projected = self.resistance_derivative_wrt_design_base * self.design_parameter_projected_derivatives
        self.resistance_derivative_wrt_design = resistance_derivative_wrt_design_projected
        # self.resistance_derivative_wrt_design[mask] = self._ApplyDiffusiveFilterDerivative(resistance_derivative_wrt_design_projected)
    
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
        self._RunStageSolutionLoop() 

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
        self.MpiPrint("--|" + top_opt_stage_str + "| START SOLUTION LOOP")
        while self.KeepAdvancingSolutionLoop():
            self.MpiPrint("--|" + top_opt_stage_str + "| ADVANCE TIME")
            self.time = self._AdvanceTime()
            self.MpiPrint("--|" + top_opt_stage_str + "| INITIALIZE SOLUTION STEP")
            self.InitializeSolutionStep()
            self.MpiPrint("--|" + top_opt_stage_str + "| UPDATE PHYSICS PARAMETERS STEP")
            self.UpdatePhysicsParametersVariablesAndSynchronize()
            self.MpiPrint("--|" + top_opt_stage_str + "| PREDICT")
            self._GetSolver().Predict()
            self.MpiPrint("--|" + top_opt_stage_str + "| SOLVE SOLUTION STEP")
            is_converged = self._GetSolver().SolveSolutionStep()
            self.MpiPrint("--|" + top_opt_stage_str + "| CHECK CONVERGENCE")
            self.__CheckIfSolveSolutionStepReturnsAValue(is_converged)
            self.MpiPrint("--|" + top_opt_stage_str + "| FINALIZE SOLUTION STEP")
            self.FinalizeSolutionStep()
        self.MpiPrint("--|" + top_opt_stage_str + "| END SOLUTION LOOP")

    def KeepAdvancingSolutionLoop(self):
        """This method specifies the stopping criteria for breaking the solution loop
        It can be overridden by derived classes
        """
        if self.IsPhysicsStage(): # NS
            return self.time < self.end_time
        elif self.IsAdjointStage(): #ADJ
            return self.time > self.start_time
        else:
            self.MpiPrint("--|" + self.topology_optimization_stage_str + "| Time check outside NS or ADJ_NS solution")
            return False

    def _AdvanceTime(self):
        """
        This methods hadls the time for the Topology Optimization problems, currently its purpos it's just to make thinks work
        """
        time = super()._AdvanceTime()
        if self.IsAdjointStage():
            time = 0.0
        return time

    def _ReInitializePhysics(self):
        """
        This method Re-Initializes the physical values in order to be consistent with the solution of the current optimization step
        """
        self.MpiPrint("--|" + self.topology_optimization_stage_str + "| RE-INITIALIZE PHYSICS")
        self._GetComputingModelPart().ProcessInfo.SetValue(KratosMultiphysics.TIME, 0.0)

    def _EvaluateOptimizationProblem(self, design_parameter = [], print_results = False):
        self.MpiPrint("\n--|EVALUATE OPTIMIZATION PROBLEM|")
        if (len(design_parameter) != 0):
            self._UpdateDesignParameterAndPhysicsParameters(design_parameter)
            self._EvaluateRequiredGradients()
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
        This method computes the resistance functional: int_{Omega}{\\alpha||u||^2}
        """
        mp = self._GetComputingModelPart()
        velocity = np.asarray(KratosMultiphysics.VariableUtils().GetSolutionStepValuesVector(mp.GetCommunicator().LocalMesh().Nodes, KratosMultiphysics.VELOCITY, 0, self.dim)).reshape(self.n_nodes, self.dim)
        nodal_velocity_norm = np.linalg.norm(velocity, axis=1)
        integrand = self.resistance * (nodal_velocity_norm**2) #component-wise multiplication
        self.functionals[0] = np.dot(self.nodal_domain_sizes, integrand)
        if (self.first_iteration):
            self.initial_functionals_values[0] = self.functionals[0] 
        if (print_functional):
            self.MpiPrint("--|" + self.topology_optimization_stage_str + "| ---> Resistance Functional (no weight):", self.functionals[0])
        else:
            self.MpiPrint("--|" + self.topology_optimization_stage_str + "| ---> Resistance Functional")

    def _EvaluateStrainRateFunctional(self, print_functional=False):
        """
        This method computes the Strain-Rate functional: int_{Omega}{\\2*mu*||1/2*[grad(u)+grad(u)^T]||^2}
        """
        mp = self._GetComputingModelPart()
        velocity = np.asarray(KratosMultiphysics.VariableUtils().GetSolutionStepValuesVector(mp.GetCommunicator().LocalMesh().Nodes, KratosMultiphysics.VELOCITY, 0, self.dim)).reshape(self.n_nodes, self.dim)
        vel = velocity[self.element_nodes_ids[:]]
        vel_gradient = np.matmul(np.transpose(vel, axes=(0,2,1)), self.shape_functions_derivatives)
        vel_symmetric_gradient = 1.0/2.0 * (vel_gradient+(np.transpose(vel_gradient, axes=(0,2,1))))
        vel_symmetric_gradient_norm_squared = (np.linalg.norm(vel_symmetric_gradient, ord='fro', axis=(1, 2)))**2
        mu = self._GetComputingModelPart().Elements[1].Properties.GetValue(KratosMultiphysics.DYNAMIC_VISCOSITY)
        self.functionals[1] = 2*mu* np.dot(vel_symmetric_gradient_norm_squared, self.elemental_domain_size)
        if (self.first_iteration):
            self.initial_functionals_values[1] = self.functionals[1]
        if (print_functional):
            self.MpiPrint("--|" + self.topology_optimization_stage_str + "| ---> Strain-Rate Functional (no weight):", self.functionals[1])
        else:
            self.MpiPrint("--|" + self.topology_optimization_stage_str + "| ---> Strain-Rate Functional")

    def _EvaluateVorticityFunctional(self, print_functional=False):
        """
        This method computes the Vorticity functional: int_{Omega}{\\2*mu*||1/2*[grad(u)-grad(u)^T]||^2}
        """
        mp = self._GetComputingModelPart()
        velocity = np.asarray(KratosMultiphysics.VariableUtils().GetSolutionStepValuesVector(mp.GetCommunicator().LocalMesh().Nodes, KratosMultiphysics.VELOCITY, 0, self.dim)).reshape(self.n_nodes, self.dim)
        vel = velocity[self.element_nodes_ids[:]]
        vel_gradient = np.matmul(np.transpose(vel, axes=(0,2,1)), self.shape_functions_derivatives)
        vel_antisymmetric_gradient = 1.0/2.0 * (vel_gradient-(np.transpose(vel_gradient, axes=(0,2,1))))
        vel_antisymmetric_gradient_norm_squared = (np.linalg.norm(vel_antisymmetric_gradient, ord='fro', axis=(1, 2)))**2
        mu = self._GetComputingModelPart().Elements[1].Properties.GetValue(KratosMultiphysics.DYNAMIC_VISCOSITY)
        self.functionals[2] = 2*mu* np.dot(vel_antisymmetric_gradient_norm_squared, self.elemental_domain_size)
        if (self.first_iteration):
            self.initial_functionals_values[2] = self.functionals[2]
        if (print_functional):
            self.MpiPrint("--|" + self.topology_optimization_stage_str + "| ---> Vorticity Functional: (no weight)", self.functionals[2])    
        else:
            self.MpiPrint("--|" + self.topology_optimization_stage_str + "| ---> Vorticity Functional")
    
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
        mp = self._GetComputingModelPart()
        velocity = np.asarray(KratosMultiphysics.VariableUtils().GetSolutionStepValuesVector(mp.GetCommunicator().LocalMesh().Nodes, KratosMultiphysics.VELOCITY, 0, self.dim)).reshape(self.n_nodes, self.dim)
        return self.functional_weights[0]*self.resistance_derivative_wrt_design_base * np.sum(velocity*velocity, axis=1) * self.nodal_domain_sizes 

    def _ComputeFunctionalDerivativesFluidPhysicsContribution(self):
        mp = self._GetComputingModelPart()
        velocity = np.asarray(KratosMultiphysics.VariableUtils().GetSolutionStepValuesVector(mp.GetCommunicator().LocalMesh().Nodes, KratosMultiphysics.VELOCITY, 0, self.dim)).reshape(self.n_nodes, self.dim)
        velocity_adjoint= np.asarray(KratosMultiphysics.VariableUtils().GetSolutionStepValuesVector(mp.GetCommunicator().LocalMesh().Nodes, KratosMultiphysics.VELOCITY_ADJ, 0, self.dim)).reshape(self.n_nodes, self.dim)        
        return self.resistance_derivative_wrt_design_base * np.sum(velocity*velocity_adjoint, axis=1) * self.nodal_domain_sizes

    def _EvaluateWSSConstraintAndDerivative(self):
        """
        For simplicity of notation we refer to the design parameter gradient as g.
        g: design parameter gradient
        gn: design parameter gradient norm 
        gn_int: integral over the domain of gn, used to normalize gn
        w: gn / gn_int (integral weights based on design parameter gradient)
        """
        mp = self._GetComputingModelPart()
        g = np.asarray(KratosMultiphysics.VariableUtils().GetSolutionStepValuesVector(mp.GetCommunicator().LocalMesh().Nodes, KratosMultiphysics.DESIGN_PARAMETER_GRADIENT, 0, self.dim)).reshape(self.n_nodes, self.dim)
        gn = np.linalg.norm(g , axis=1)
        gn_max = np.max(gn)
        if (gn_max > 1e-14):
            gn_int = np.dot(gn, self.nodal_domain_sizes)
            v_grad = self._AssembleVelocityGradientOnNodes()
            nodal_tangents_to_g = self._EvaluateOrthogonalBasis(g)
            # prod = np.einsum('ijk,ik->ij',  nodal_tangents_to_g, g)
            psi_vect =np.einsum('ijk,ik,ilj->il', v_grad, g, nodal_tangents_to_g)
            psi =  np.linalg.norm(psi_vect, axis=1)
            mu = self._GetComputingModelPart().Elements[1].Properties.GetValue(KratosMultiphysics.DYNAMIC_VISCOSITY)
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
        self.MpiPrint("\n--|PRINT OPTIMIZATION PROBLEM DATA|")
        self._PrintFunctionals()
        self._PrintConstraints()
    
    def _PrintFunctionals(self):
        self.MpiPrint("--|" + self.topology_optimization_stage_str + "| TOTAL FUNCTIONAL  : " +  str(self.functional))
        self.MpiPrint("--|" + self.topology_optimization_stage_str + "| INITIAL FUNCTIONAL: " +  str(self.initial_functional))
        if (abs(self.functional_weights[0]) > 1e-10):
            self.MpiPrint("--|" + self.topology_optimization_stage_str + "| ---> Resistance Functional (" + str(self.functional_weights[0]) + "): " + str(self.weighted_functionals[0]/self.initial_functional_abs_value))
        if (abs(self.functional_weights[1]) > 1e-10):
            self.MpiPrint("--|" + self.topology_optimization_stage_str + "| ---> Strain-Rate Functional (" + str(self.functional_weights[1]) + "): " + str(self.weighted_functionals[1]/self.initial_functional_abs_value))
        if (abs(self.functional_weights[2]) > 1e-10):
            self.MpiPrint("--|" + self.topology_optimization_stage_str + "| ---> Vorticity Functional (" + str(self.functional_weights[2]) + "): " + str(self.weighted_functionals[2]/self.initial_functional_abs_value))
        
    def _PrintFunctionalsToFile(self):
        with open("functional_history.txt", "a") as file:
            file.write(str(self.functional) + " ")
            for ifunc in range(self.n_functionals): 
                file.write(str(self.weighted_functionals[ifunc]) + " ")
            file.write("\n")

    def _EvaluateDesignParameterChange(self):
        old_design_parameter_norm = np.linalg.norm(1-self.old_design_parameter)
        self.design_parameter_change = np.linalg.norm(self.design_parameter-self.old_design_parameter)
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
        self.MpiPrint("--|" + self.topology_optimization_stage_str + "| DESIGN PARAMETER CHANGE: " + str(self.design_parameter_change))
        return design_parameter_converged        

    def _ResetFunctionalOutput(self):
        with open("functional_history.txt", "w") as file:
            file.write("TOTAL | RESISTANCE | STRAIN-RATE | VORTICITY | OUTLET CONCENTRATION | NONE\n")
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
        self.MpiPrint("--|" + self.topology_optimization_stage_str + "| VOLUME FRACTION: " + str(self.volume_fraction))
        self.MpiPrint("--|" + self.topology_optimization_stage_str + "| ---> Volume Constraint: " + str(self.volume_constraint))

    def _PrintWSSConstraint(self):
        self.MpiPrint("--|" + self.topology_optimization_stage_str + "| WSS VALUE: " + str(self.wss_value))
        self.MpiPrint("--|" + self.topology_optimization_stage_str + "| ---> WSS Resistance: " + str(self.resistance_parameters["value_full"].GetDouble()))
        self.MpiPrint("--|" + self.topology_optimization_stage_str + "| ---> WSS Constraint: " + str(self.wss_constraint))

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
        if (self.apply_diffusive_filter):
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
        only_opt_mp_nodes = self._GetModelPartNodesSubset(mp, self.nodes_ids_local_partition_to_global_dictionary[mask])
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
        self.pde_diffusive_filter_process = ApplyTopologyOptimizationPdeFilterProcessMpi(self.model, self.diffusive_filter_type_settings, self.diffusive_filter_radius, self._GetMainModelPart(),  self._GetOptimizationDomain(), self._GetOptimizationDomainNodesMask(), self.nodes_ids_global_to_local_partition_dictionary)

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
        if (self.apply_diffusive_filter):
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
        if (self.apply_diffusive_filter):
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
        mp = self._GetComputingModelPart()
        design_parameter = np.zeros(self.n_nodes)
        count = 0
        for node in mp.GetCommunicator().LocalMesh().Nodes:
            design_parameter[count] = node.GetSolutionStepValue(KratosMultiphysics.DESIGN_PARAMETER)
            count += 1
        return design_parameter
    
    def OutputSolutionStep(self):
        self.MpiPrint("\n--| PRINT SOLUTION STEP OUTPUT TO FILES")
        super().OutputSolutionStep()

    def _PrintSolution(self):
        self.OutputSolutionStep()
        self._PrintFunctionalsToFile()

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
                self.local_hessian = KratosMMG.ComputeHessianSolMetricProcess2D(main_mp, KratosMultiphysics.DISTANCE, metric_param)
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
        if (abs(self.functional_weights[0]) > 1e-10):
            self._EvaluateResistanceFunctional(print_functional)
        if (abs(self.functional_weights[1]) > 1e-10):
            self._EvaluateStrainRateFunctional(print_functional)
        if (abs(self.functional_weights[2]) > 1e-10):
            self._EvaluateVorticityFunctional(print_functional)

    def _CheckMaterialProperties(self, check = False):
        if (check):
            self.MpiPrint("--|CHECK| Check Fluid Properties")
            self._GetSolver()._CheckMaterialProperties()

    def _UpdateRelevantPhysicsVariables(self):
        pass

    def _UpdateRelevantAdjointVariables(self):
        pass

    def _EvaluateRequiredGradients(self):
        self._ComputeScalarVariableNodalGradient(KratosMultiphysics.DESIGN_PARAMETER, KratosMultiphysics.DESIGN_PARAMETER_GRADIENT)
        if (self.use_wss_constraint):
            self._ComputeScalarVariableNodalGradient(KratosMultiphysics.VELOCITY_X, KratosMultiphysics.VELOCITY_X_GRADIENT)
            self._ComputeScalarVariableNodalGradient(KratosMultiphysics.VELOCITY_Y, KratosMultiphysics.VELOCITY_Y_GRADIENT)
            if (self.dim == 3):
                self._ComputeScalarVariableNodalGradient(KratosMultiphysics.VELOCITY_Z, KratosMultiphysics.VELOCITY_Z_GRADIENT)

    def _AssembleVelocityGradientOnNodes(self):
        mp = self._GetComputingModelPart()
        gradient_x = np.asarray(KratosMultiphysics.VariableUtils().GetSolutionStepValuesVector(mp.GetCommunicator().LocalMesh().Nodes, KratosMultiphysics.VELOCITY_X_GRADIENT, 0, self.dim)).reshape(self.n_nodes, self.dim)
        gradient_y = np.asarray(KratosMultiphysics.VariableUtils().GetSolutionStepValuesVector(mp.GetCommunicator().LocalMesh().Nodes, KratosMultiphysics.VELOCITY_Y_GRADIENT, 0, self.dim)).reshape(self.n_nodes, self.dim)
        gradient = np.zeros((self.n_nodes, self.dim, self.dim))
        gradient[:,0,:] = gradient_x
        gradient[:,1,:] = gradient_y
        if (self.dim == 3):
            gradient_z = np.asarray(KratosMultiphysics.VariableUtils().GetSolutionStepValuesVector(mp.GetCommunicator().LocalMesh().Nodes, KratosMultiphysics.VELOCITY_Z_GRADIENT, 0, self.dim)).reshape(self.n_nodes, self.dim)
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
            "resistance": {
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
            }
        }""")
        default_physics_parameters_settings.AddMissingParameters(self.GetBasePhysicsParametersSettings())
        return default_physics_parameters_settings
    
    def GetBasePhysicsParametersSettings(self):
        ##settings string in json format
        base_physics_parameters_settings = KratosMultiphysics.Parameters("""
        {
            "resistance": {},                                                                
            "conductivity": {},
            "decay": {},
            "convection_coefficient": {},
            "transport_source": {}                                                     
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
                }
            },
            "optimization_problem_settings": {
                "functional_weights": {
                    "fluid_functionals": {
                        "resistance" : {
                            "weight": 1.0
                        },
                        "strain_rate" : {
                            "weight": 1.0
                        },
                        "vorticity" : {
                            "weight": 0.0
                        }
                    },
                    "transport_functionals":
                    {
                        "outlet_transport_scalar" : {
                            "weight"    : 0.0,
                            "outlet_model_part": "Outlet",
                            "target_value": 0.0
                        },
                        "focus_region_transport_scalar" : {
                            "weight"    : 0.0,
                            "focus_region_model_part": "FocusRegion",
                            "target_value": 0.0
                        },
                        "conductivity_transfer" : {
                            "weight": 1.0
                        },
                        "convection_transfer" : {
                            "weight": 0.0
                        },
                        "decay_transfer" : {
                            "weight": 0.0
                        },
                        "source_transfer" : {
                            "weight": 0.0
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
                },
                "use_custom_initial_design": false,
                "custom_initial_design_settings": { 
                    "custom_initial_design": [{
                        "model_part_name": "",
                        "initial_value": "0.0"
                    }]
                }
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
            "remeshing_settings": {},
            "symmetry_settings": {},
            "optimization_domain_settings": {},
            "optimization_problem_settings": {}                                                  
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
        super().__init__(model,parameters) 
        self._ReadOptimizationParameters()
        # self._CreateTopologyOptimizationSolvers() # currently it is a useless method 
        self._SetMinMaxIt()  
        self._SetTopologyOptimizationName()
        self.InitializeDataCommunicator()

    def InitializeDataCommunicator(self):
        self.data_communicator = DataCommunicator.GetDefault()

    def PrepareAdjointSolver(self):
        """This method prepares the Adjoint Navier-Stokes Solver in the AnalysisStage 
        Usage: It is designed to be called ONCE, BEFORE the execution of the solution-loop
        Prepare Solver : ImportModelPart -> PrepareModelPart -> AddDofs
        """
        self._GetAdjointSolver().ImportModelPart(model_parts=self._GetPhysicsMainModelPartsList(), physics_solver_distributed_model_part_importer=self.physics_solver.distributed_model_part_importer)
        self._GetAdjointSolver().PrepareModelPart()
        self._GetAdjointSolver().AddDofs()

    def _CreateNodesIdsDictionary(self):
        self.MpiPrint("--|" + self.topology_optimization_stage_str + "| ---> Creates Nodes Ids Dictionary")
        self.nodes_ids_global_to_local_partition_dictionary = {}
        self.nodes_ids_local_partition_to_global_dictionary = {}
        count = 0
        mp = self._GetComputingModelPart()
        for node in mp.GetCommunicator().LocalMesh().Nodes:
            self.nodes_ids_global_to_local_partition_dictionary[node.Id] = count
            self.nodes_ids_local_partition_to_global_dictionary[count]   = node.Id
            count += 1
        if (count != len(self.nodes_ids_global_to_local_partition_dictionary.keys())):
            raise RuntimeError("Wrong reordering of nodes ids. The counted number of nodes is different from len(mp.GetCommunicator().LocalMesh().Nodes).")
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
        # 0: resistance  : int_{Omega}{alpha*||u||^2}
        # 1: strain-rate : int_{Omega}{2*mu*||S||^2} , with S = 1/2*(grad(u)+grad(u)^T) strain-rate tensor
        # 2: vorticity   : int_{Omega}{2*mu*||R||^2} = int_{Omega}{mu*||curl(u)||^2} , curl(u) = vorticity vector, R = 1/2*(grad(u)-grad(u)^T) rotation-rate tensor
        # 3: outlet_transport_scalar : int_{Gamma_{out}}{c}
        # 4: region_transport_scalar: int_{Omega}{c^2}
        # 5: transport_scalar_diffusion: int_{Omega}{D\\||grad(u)||^2}
        # 6: transport_scalar_convection: int_{Omega}{beta*T*dot(u,grad(T))}
        # 7: transport_scalar_decay: int_{Omega}{kT^2}
        # 8: transport_scalar_source: int_{Omega}{-Q*T}
        """
        self._SetTopologyOptimizationStage(3)
        self.MpiPrint("--|" + self.topology_optimization_stage_str + "| EVALUATE FUNCTIONAL VALUE")
        self.EvaluateFunctionals(print_functional)
        self.functional = np.dot(self.functional_weights, self.functionals)
        self.weighted_functionals = self.functional_weights * self.functionals
        if _CheckIsDistributed():
            self.MpiSynchronizeLocalFunctionalValues()
        if (self.MpiRunOnlyRank(0)):
            if (self.first_iteration):
                    self.initial_functional = self.functional
                    self.initial_functional_abs_value = abs(self.initial_functional)
                    if (abs(self.initial_functional_abs_value) < 1e-10):
                        self.initial_functional_value = 1.0
                        self.initial_functional_abs_value = 1.0
                    self.initial_functionals_abs_value = np.abs(self.functionals)
                    self.initial_weighted_functionals_abs_value = np.abs(self.weighted_functionals)
                    # self.first_iteration = False
            self.functional = self.functional / self.initial_functional_abs_value

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
        if _CheckIsDistributed():
            # self.MpiBarrier()
            local_design_parameter_integral = np.dot(self.design_parameter, self.nodal_optimization_domain_sizes)
            total_design_parameter_integral = self.data_communicator.SumAll(local_design_parameter_integral)
            self.design_parameter_integral  = total_design_parameter_integral
        else:
            self.design_parameter_integral  = np.dot(self.design_parameter, self.nodal_optimization_domain_sizes)

    def _PrintOptimizationProblem(self):
            # self.MpiBarrier()
            if (self.MpiRunOnlyRank(0)):
                self.MpiPrint("\n--|PRINT OPTIMIZATION PROBLEM DATA|")
                self._PrintFunctionals()
                self._PrintConstraints()

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
        self.nodal_optimization_domain_sizes = np.zeros(self.n_opt_design_parameters)
        mp = self._GetOptimizationDomain()
        nodal_area_process = KratosMultiphysics.CalculateNodalAreaProcess(mp, self.dim)
        nodal_area_process.Execute()
        self._UpdateNodalOptimizationDomainSizeArrayFromNodalAreaVariable()
        self._CorrectNodalOptimizationDomainSizeWithSymmetry()

    def _UpdateNodalDomainSizeArrayFromNodalAreaVariable(self):
        mp = self._GetComputingModelPart()
        for node in mp.GetCommunicator().LocalMesh().Nodes:
            current_rank_node_id = self.nodes_ids_global_to_local_partition_dictionary[node.Id]
            self.nodal_domain_sizes[current_rank_node_id] = node.GetSolutionStepValue(KratosMultiphysics.NODAL_AREA)

    def _UpdateNodalOptimizationDomainSizeArrayFromNodalAreaVariable(self):
        mp = self._GetOptimizationDomain()
        for node in mp.GetCommunicator().LocalMesh().Nodes:
            current_rank_optimization_node_id = self.global_to_opt_mask_nodes[self.nodes_ids_global_to_local_partition_dictionary[node.Id]]
            self.nodal_optimization_domain_sizes[current_rank_optimization_node_id] = node.GetSolutionStepValue(KratosMultiphysics.NODAL_AREA)

    def _ComputeOptimizationDomainNodesMask(self):
        """
        This method build the mmask to pass from the optimization domain nodes to the total domain nodes
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
            for node in non_opt_mp.GetCommunicator().LocalMesh().Nodes:
                self.global_to_opt_mask_nodes[self.nodes_ids_global_to_local_partition_dictionary[node.Id]] = -1
            n_non_opt_nodes = len(non_opt_mp.GetCommunicator().LocalMesh().Nodes)
        else:
            n_non_opt_nodes = 0
        self.n_opt_design_parameters = self.n_nodes - n_non_opt_nodes
        self.optimization_domain_nodes_mask = np.zeros(self.n_opt_design_parameters, dtype=int)
        count_opt_nodes = 0
        for node in mp.GetCommunicator().LocalMesh().Nodes:
            full_domain_to_local_rank_node_id = self.nodes_ids_global_to_local_partition_dictionary[node.Id]
            if (self.global_to_opt_mask_nodes[full_domain_to_local_rank_node_id] != -1): #if the node is only in the optimization domain, add it to the mask
                self.optimization_domain_nodes_mask[count_opt_nodes] = full_domain_to_local_rank_node_id
                self.global_to_opt_mask_nodes[full_domain_to_local_rank_node_id] = count_opt_nodes
                count_opt_nodes +=1
        if (count_opt_nodes != self.n_opt_design_parameters):
            self.MpiPrint("!!! WARNING: wrong initialization of the Optimization Domain Nodes Mask")

    def _UpdateFunctionalDerivativesVariable(self):
        mp = self._GetComputingModelPart()
        for node in mp.GetCommunicator().LocalMesh().Nodes:
        #     node.SetValue(KratosMultiphysics.FUNCTIONAL_DERIVATIVE, self.functional_derivatives_wrt_design[node.Id-1][0])
        # # synchronize the value across all the ranks of the data_communicator
        # # self.MpiBarrier()
        # mp = self._GetComputingModelPart()
        # for node in mp.GetCommunicator().LocalMesh().Nodes:
        #     if node.Is(KratosMultiphysics.INTERFACE):  # Only shared nodes
        #         local_value = node.GetValue(KratosMultiphysics.FUNCTIONAL_DERIVATIVE)
        #         summed_value = self.data_communicator.SumAll(local_value)  # Sum across all ranks
        #         node.SetValue(KratosMultiphysics.FUNCTIONAL_DERIVATIVE, summed_value)  # Store the summed result
            node.SetValue(KratosMultiphysics.FUNCTIONAL_DERIVATIVE, self.functional_derivatives_wrt_design[self.nodes_ids_global_to_local_partition_dictionary[node.Id]][0])
        
    def MpiSynchronizeLocalFunctionalValues(self):
        # self.MpiBarrier()
        local_values = self.functionals
        # Sum the values across all ranks
        total_values = self.data_communicator.SumAll(local_values)
        self.functionals = total_values
        if (self.MpiRunOnlyRank(0)):
            self.weighted_functionals  = self.functional_weights * self.functionals
            self.functional  = np.dot(self.functional_weights, self.functionals)
        # self.MpiBarrier()

    def UpdatePhysicsParametersVariablesAndSynchronize(self):
        self.UpdatePhysicsParametersVariables()
        self._SynchronizePhysicsParametersVariables()
        
    def _SynchronizePhysicsParametersVariables(self):
        self._GetMainModelPart().GetCommunicator().SynchronizeVariable(KratosMultiphysics.RESISTANCE)

    def _SolveMMA(self, design_parameter, n_opt_variables, n_opt_constraints, min_value, max_value, max_outer_it, kkt_tolerance):
        # self.MpiBarrier()
        self.MpiPrint("--|" + self.topology_optimization_stage_str + "| SOLVE MMA")
        # MMA PARAMETERS INITIALIZATION
        # self.MpiBarrier()
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
        local_value = local_value.flatten().tolist()
        gather_local_value = self.data_communicator.AllGathervDoubles(local_value)
        for irank in range(n_ranks):
            gather_local_value[irank] = np.asarray(gather_local_value[irank]).reshape(n_opt_variables_in_rank[irank],1)
        result = np.concatenate(gather_local_value, axis=0)
        return result
    
    def CreateConstraintsDerivativesCompleteArrays(self, local_value, n_ranks, n_opt_variables_in_rank, n_opt_constraints):
        local_value = local_value.flatten().tolist()
        gather_local_value = self.data_communicator.AllGathervDoubles(local_value)
        for irank in range(n_ranks):
            gather_local_value[irank] = np.asarray(gather_local_value[irank]).reshape(n_opt_constraints, n_opt_variables_in_rank[irank])
        result = np.concatenate(gather_local_value, axis=1)
        return result
    
    def CreateFunctionalAndConstraintsDerivativesCompleteArrays(self, temp_df0dx, temp_dfdx, n_ranks, n_opt_variables_in_rank, n_opt_constraints):
        # self.MpiBarrier()
        df0dx = self.CreateFunctionalDerivativesCompleteArrays(temp_df0dx, n_ranks, n_opt_variables_in_rank)
        # self.MpiBarrier()
        dfdx = self.CreateConstraintsDerivativesCompleteArrays(temp_dfdx, n_ranks, n_opt_variables_in_rank, n_opt_constraints)
        # self.MpiBarrier()
        return df0dx, dfdx

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
        if (not _CheckIsDistributed()):
            return True
        elif (self.data_communicator.Rank() == rank):
            return True
        else:
            return False
        
    def MpiPrint(self, text_to_print="", rank=0, set_barrier=False):
        if (not _CheckIsDistributed()):
            print(text_to_print)
        else:
            if (set_barrier):
                self.MpiBarrier()
            if (self.MpiRunOnlyRank(rank)):
                print(text_to_print)
            if (set_barrier):
                self.MpiBarrier()    
