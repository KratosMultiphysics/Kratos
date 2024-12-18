from sys import argv

import numpy as np #import the numpy library
import scipy as sp #import the scipy library
from scipy.spatial import KDTree
from scipy.sparse import dok_matrix, lil_matrix
import mmapy as MMA # import the MMA subroutines python library: https://github.com/arjendeetman/GCMMA-MMA-Python
import time

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
        self._SetMinMaxIt()  

    def _CreateTopologyOptimizationSolvers(self):
        """
        This method creates the NS and ADJ_NS solvers
        """
        self.NS_solver = fluid_topology_optimization_solver.CreateSolver(self.model, self.project_parameters)
        self.ADJ_NS_solver = fluid_topology_optimization_solver.CreateSolver(self.model, self.project_parameters, True)

    def _CreateSolver(self, isAdjointSolver = False):
        """
        This method creates a solver
        isAdjointSolver == False --> NS_solver
        isAdjointSolver == True  --> ADJ_NS_solver
        """
        return fluid_topology_optimization_solver.CreateSolver(self.model, self.project_parameters, isAdjointSolver)
    
    def _GetSolver(self, force_adjoint = False):
        """
        This method returns solver in use for the current topology optimization phase
        """
        if not hasattr(self, 'NS_solver'):
                self.NS_solver = self._CreateSolver()
        if not hasattr(self, 'ADJ_NS_solver'):
                self.ADJ_NS_solver = self._CreateSolver(True)
        self._solver = self._GetTopologyOptimizationStageSolver(force_adjoint)
        return self._solver
    
    def _GetNavierStokesSolver(self):
        """
        This method returns the Navier-Stokes Solver
        """
        if not hasattr(self, 'NS_solver'):
                self.NS_solver = self._CreateSolver(False)
        return self.NS_solver
    
    def _GetAdjointNavierStokesSolver(self):
        """
        This method returns the Adjoint Navier-Stokes Solver
        """
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
        
    def _SetMinMaxIt(self, min_iteration = 0, max_iteration = 1):
        """
        This method sets the minimum & maximum iteration for the topology optimization process
        """
        self.min_it = min_iteration
        self.max_it = max_iteration
        if (self.min_it > self.max_it):
            print("\n!!!WARNING: wrong initialization of the min & max number of iterations\n")

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
        This method executes the entire Fluid Topology Optimization Stage
        """
        self.Initialize()
        self.RunSolutionLoop()
        self.Finalize()

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
        if self._GetComputingModelPart().ProcessInfo[Kratos.IS_RESTARTED]:
            self.time = self._GetComputingModelPart().ProcessInfo[Kratos.TIME]
        else:
            self.time = self.project_parameters["problem_data"]["start_time"].GetDouble()
            self._GetComputingModelPart().ProcessInfo[Kratos.TIME] = self.time
        self.start_time = self.time
        self.end_time = self.project_parameters["problem_data"]["end_time"].GetDouble()
        ## If the echo level is high enough, print the complete list of settings used to run the simulation
        if self.echo_level > 1:
            with open("ProjectParametersOutput.json", 'w') as parameter_output_file:
                parameter_output_file.write(self.project_parameters.PrettyPrintJsonString())
        # Print Start Analysis
        Kratos.Logger.PrintInfo("\n" + self._GetSimulationName(), "Analysis -START- ")

    def PrepareSolvers(self):
        """This method prepares the NS and ADJ_NS Solvers in the AnalysisStage 
        Usage: It is designed to be called ONCE, BEFORE the execution of the solution-loop
        Prepare Solvers : ImportModelPart -> PrepareModelPart -> AddDofs
        """
        self.PrepareNavierStokesSolver()
        self.PrepareAdjointNavierStokesSolver()
    
    def PrepareNavierStokesSolver(self):
        """This method prepares the Navier-Stokes primal problem Solver in the AnalysisStage 
        Usage: It is designed to be called ONCE, BEFORE the execution of the solution-loop
        Prepare Solver : ImportModelPart -> PrepareModelPart -> AddDofs
        """
        self._GetNavierStokesSolver().ImportModelPart()
        self._GetNavierStokesSolver().PrepareModelPart()
        self._GetNavierStokesSolver().AddDofs()
    
    def PrepareAdjointNavierStokesSolver(self):
        """This method prepares the Adjoint Navier-Stokes Solver in the AnalysisStage 
        Usage: It is designed to be called ONCE, BEFORE the execution of the solution-loop
        Prepare Solver : ImportModelPart -> PrepareModelPart -> AddDofs
        """
        # Modelers:
        self._GetAdjointNavierStokesSolver().ImportModelPart()
        self._GetAdjointNavierStokesSolver().PrepareModelPart()
        self._GetAdjointNavierStokesSolver().AddDofs()

    def _SetFunctionalWeights(self, weights = [1, 1, 0, 0, 0]):
        self.n_functionals = len(weights)
        self.initial_functionals_values = np.zeros(self.n_functionals)
        self.functionals = np.zeros(self.n_functionals)
        self.functional_weights = self._NormalizeFunctionalWeights(np.asarray(weights))
        self._GetComputingModelPart().ProcessInfo.SetValue(KratosCFD.FUNCTIONAL_WEIGHTS, self.functional_weights)

    def _NormalizeFunctionalWeights(self, weights):
        weights_sum = np.sum(np.abs(weights))
        return np.asarray(weights)/weights_sum
    
    def InitializeSolvers(self):
        """This method initializes the NS and ADJ_NS Solvers in the AnalysisStage 
        Usage: It is designed to be called ONCE, BEFORE the execution of the solution-loop
        """
        self.InitializeNavierStokesSolver()
        self.InitializeAdjointNavierStokesSolver()

    def InitializeNavierStokesSolver(self):
        """This method initializes the NS Solver in the AnalysisStage 
        Usage: It is designed to be called ONCE, BEFORE the execution of the solution-loop
        """
        self._GetNavierStokesSolver().Initialize()

    def InitializeAdjointNavierStokesSolver(self):
        """This method initializes the ADJ_NS Solver in the AnalysisStage 
        Usage: It is designed to be called ONCE, BEFORE the execution of the solution-loop
        """
        self._GetAdjointNavierStokesSolver().Initialize()

    def Check(self):
        """This method checks the AnalysisStage
        Keyword arguments:
        self -- It signifies an instance of a class.
        """
        # Checking solver
        self._GetNavierStokesSolver().Check()
        self._GetAdjointNavierStokesSolver().Check()

        # Checking processes
        for process in self._GetListOfProcesses():
            process.Check()

    def RunSolutionLoop(self):
        self._InitializeTopologyOptimizationProblem()
        self.SolveFluidTopologyOptimization()

    def _InitializeTopologyOptimizationProblem(self):
        """
        This method Initializes the topology optimization problem solution
        """
        print("\n--------------------------------------------------------")
        print(  "--| FLUID TOPOLOGY OPTIMIZATION PREPROCESSING")
        print("--------------------------------------------------------")
        print("--|INITIALIZE|")
        self._ResetFunctionalOutput()
        self._GeometricalPreprocessing()
        self._InitializeOptimization()
        self._InitializeDomainDesign()
        self._PreprocessDerivatives() 

    def SolveFluidTopologyOptimization(self):
        self.design_parameter_change = self.design_parameter_change_toll + 10
        end_solution = False
        while (not end_solution):
            self.opt_it = self.opt_it+1
            print("\n--------------------------------------------------------")
            print(  "--| FLUID TOPOLOGY OPTIMIZATION SOLUTION LOOP. IT:", self.opt_it)
            print("--------------------------------------------------------")
            self.old_design_parameter = self.design_parameter
            self._SolveOptimizer(self.opt_it == 1)
            self._SolveTopologyOptimizationStepPhysics()
            self._EvaluateOptimizationProblem(first_iteration=(self.opt_it==1), print_results=True)
            self._PrintFunctionalsToFile()
            end_solution = self._IsTopologyOptimizationSolutionEnd()
            if (end_solution):
                print("\n--------------------------------------------------------")
                print("--| ENDING FLUID TOPOLOGY OPTIMIZATION SOLUTION LOOP")
                if (self.converged):
                    print("--| ---> CONVERGED!")
                else:
                    print("--| ---> Reached Max Number of Iterations")
                print("--------------------------------------------------------\n")
            
    def _IsTopologyOptimizationSolutionEnd(self):
        design_parameter_converged = self._EvaluateDesignParameterChange()
        volume_constraint_valid = not (self.volume_constraint > 0.0)
        self.converged = design_parameter_converged and volume_constraint_valid
        return not (((self.opt_it < self.max_it) and (not self.converged)) or (self.opt_it < self.min_it))

    def _SolveOptimizer(self, first_iteration = False):
        print("\n--|OPTIMIZER|")
        self._SetTopologyOptimizationStage(3)
        if (first_iteration):
            print("--|" + self.topology_optimization_stage_str + "| First Iteration: DO NOTHING")
            pass
        else:
            print("--|" + self.topology_optimization_stage_str + "| SOLVE OPTIMIZER")
            opt_design_parameter = self._ExtractVariableInOptimizationDomain(self.design_parameter)
            self._SolveMMA(opt_design_parameter, self.n_opt_design_parameters, self.n_optimization_constraints)
    
    def _SolveTopologyOptimizationStepPhysics(self): 
        """
        This method organizs the physics solution of an optimization step
        """     
        ## Current Optimization Step Solutions 
        self._InitializeTopologyOptimizationStepPhysicsSolution((self.opt_it==1))
        self._SolveNaviersStokesProblem() # NAVIER-STOKES PROBLEM SOLUTION
        # SOLVE TRANSPORT
        # SOLVE ADJ TRANSPORT
        self._SolveAdjointNaviersStokesProblem() # ADJOINT NAVIER-STOKES PROBLEM SOLUTION   
    
    def _GeometricalPreprocessing(self):
        """
        This method preprocess all the useful values and quantities for the topology optimization solution
        """
        print("--|" + self.topology_optimization_stage_str + "| GEOMETRICAL PREPROCESSING")
        mp = self._GetComputingModelPart()
        self.dim = mp.ProcessInfo.GetValue(Kratos.DOMAIN_SIZE)
        self.nodes_in_element = self.dim+1 # only for triangles(2D) and tetrahedra(3D)
        self._OrderNodes()
        self._OrderElements()
        self._ComputeNodalDomainSizes()
    
    def _InitializeOptimization(self):  
        print("--|" + self.topology_optimization_stage_str + "| INITIALIZE OPTIMIZATION")
        self._OptimizationGeometricalPreprocessing()
        self.opt_it = 0
        self._SetDesignParameterChangeTolerance()
        self.n_optimization_constraints = 1  
        self.constraints = np.zeros(self.n_optimization_constraints)  
        self.constraints_derivatives_wrt_design = np.zeros((self.n_optimization_constraints, self.n_opt_design_parameters))

    def _SetDesignParameterChangeTolerance(self, design_parameter_change_toll = 0.0):
        self.design_parameter_change_toll = design_parameter_change_toll

    def _InitializeDomainDesign(self):
        print("--|" + self.topology_optimization_stage_str + "| INITIALIZE DOMAIN DESIGN")
        self._SetMaxDomainVolumeFraction()
        self._InitializeDomainDesignParameter()
        self._InitializeResistance()

    def _SetMaxDomainVolumeFraction(self, volume_fraction = 1.0):
        self.max_volume_fraction = volume_fraction
        self.volume_fraction = 0.0

    def _InitializeDomainDesignParameter(self, initial_value = 0.0):
        """
        This method will handle the design parameter initialization across the whole domain.
        self.design_parameter is defined in every node of the mesh.
        The optimization process will update only the value of the nodes belonging to the optimization_domain msub model part
        """
        mask = self._GetOptimizationDomainNodesMask()
        self.volume_fraction = initial_value
        self.design_parameter = np.zeros(self.n_nodes) 
        self.design_parameter[mask] = initial_value
        mp = self._GetComputingModelPart()
        for node in mp.Nodes:
            node.SetValue(KratosCFD.DESIGN_PARAMETER, self.design_parameter[node.Id-1])

    def _OrderNodes(self):
        """
        This method orders the model part nodes in increasing order starting from 1
        """
        print("--|" + self.topology_optimization_stage_str + "| ---> Order Nodes")
        count = 0
        mp = self._GetComputingModelPart()
        for node in mp.Nodes:
            count = count+1
            node.Id = count
        if (count != len(mp.Nodes)):
            print("\nERROR: wrong reordering of nodes ids. The counted number of nodes is different from len(mp.Nodes)\n")
        self.n_nodes = count
        self.n_elements = len(mp.Elements)

    def _OrderElements(self):
        """
        This method orders the model part elements in increasing order starting from 1
        """
        print("--|" + self.topology_optimization_stage_str + "| ---> Order Elements")
        count = 0
        mp = self._GetComputingModelPart()
        for el in mp.Elements:
            count = count+1
            el.Id = count
        if (count != len(mp.Elements)):
            print("\nERROR: wrong reordering of elements ids. The counted number of elements is different from len(mp.Elements)\n")
        self.n_elements = count
    
    def _ComputeNodalDomainSizes(self):
        """
        This method compute the nodal domain size - not vectorized but it is done once. WORKS ONLY FOR TRIANGULAR AND TETRAHEDRAL MESH
        """
        print("--|" + self.topology_optimization_stage_str + "| ---> Compute Domain Size")
        mp = self._GetComputingModelPart()
        contribution_factor = 1.0/self.nodes_in_element 
        self.nodal_domain_sizes = np.zeros(len(mp.Nodes))
        self.elemental_domain_size = np.zeros(len(mp.Elements))
        self.total_domain_size = 0.0
        for elem in mp.Elements: 
            geom =  elem.GetGeometry()
            elem_domain_size = geom.DomainSize()
            self.elemental_domain_size[elem.Id-1] = elem_domain_size
            self.total_domain_size = self.total_domain_size + elem_domain_size
            for node in geom:
                self.nodal_domain_sizes[node.Id-1] += elem_domain_size * contribution_factor

    def _OptimizationGeometricalPreprocessing(self):
        self._ComputeOptimizationDomainSize()
        self._ComputeOptimizationDomainNodesMask()
        self._ComputeDesignParameterDiffusiveFilterUtilities()
        self._ComputeDesignParameterProjectiveFilterUtilities()
    
    def _ComputeOptimizationDomainSize(self):
        """
        This method compute the nodal optimization domain size. WORKS ONLY FOR TRIANGULAR AND TETRAHEDRAL MESH
        The idea is that the nodal contribution is defined on each node, but it is != 0.0 only in optimization nodes.
        In this way we can evaluate the volume constraint using a simple dot product with self.design_parameter
        """
        print("--|" + self.topology_optimization_stage_str + "| ---> Compute Optimization Domain Size")
        opt_mp = self._GetOptimizationDomain()
        contribution_factor = 1.0/self.nodes_in_element 
        self.nodal_optimization_domain_sizes = np.zeros(self.n_nodes)
        self.optimization_domain_size = 0.0
        for elem in opt_mp.Elements:
            geom =  elem.GetGeometry()
            elem_domain_size = geom.DomainSize()
            self.optimization_domain_size = self.optimization_domain_size + elem_domain_size
            for node in geom:
                self.nodal_optimization_domain_sizes[node.Id-1] += elem_domain_size * contribution_factor

    def _ComputeOptimizationDomainNodesMask(self):
        """
        This method build the mmask to pass from the optimization domain nodes to the total domain nodes
        """
        print("--|" + self.topology_optimization_stage_str + "| ---> Compute Optimization Domain Nodes Mask")
        non_opt_mp= self._GetNonOptimizationDomain()
        opt_mp= self._GetOptimizationDomain()
        mp = self._GetComputingModelPart()

        #define the boolean search for the nodes belonging only to the optimization domain model part
        self.is_only_opt_domain_node = np.ones(self.n_nodes, dtype=int)
        if (not (non_opt_mp is None)):
            for node in non_opt_mp.Nodes:
                self.is_only_opt_domain_node[node.Id-1] = 0
            n_non_opt_nodes = len(non_opt_mp.Nodes)
        else:
            n_non_opt_nodes = 0

        #find the nodes belonging to the interface between the opt_mp and the non_opt_mp
        # self.opt_non_opt_interface_nodes = np.zeros(self.n_nodes, dtype=int)
        # count_nodes_at_opt_non_opt_interface = 0
        # for node in opt_mp.Nodes:
        #     if (self.is_only_opt_domain_node[node.Id] == 0):
        #         self.opt_non_opt_interface_nodes[count_nodes_at_opt_non_opt_interface] = node.Id-1
        #         count_nodes_at_opt_non_opt_interface = count_nodes_at_opt_non_opt_interface+1
        # self.opt_non_opt_interface_nodes = self.opt_non_opt_interface_nodes[:count_nodes_at_opt_non_opt_interface]
        # self.n_nodes_at_opt_non_opt_interface = len(self.opt_non_opt_interface_nodes)

        # comput the mask to bass between the global domain nodes and the optimization_domain ones
        self.n_opt_design_parameters = self.n_nodes - n_non_opt_nodes
        self.optimization_domain_nodes_mask = np.zeros(self.n_opt_design_parameters, dtype=int)
        count_opt_nodes = 0
        for node in mp.Nodes:
            if (self.is_only_opt_domain_node[node.Id-1] == 1): #if the node is only in the optimization domain, add it to the mask
                self.optimization_domain_nodes_mask[count_opt_nodes] = node.Id-1
                count_opt_nodes = count_opt_nodes+1
        if (count_opt_nodes != self.n_opt_design_parameters):
            print("!!! WARNING: wrong initialization of the Optimization Domain Nodes Mask")

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
        print("--|" + self.topology_optimization_stage_str + "| PREPROCESS DERIVATIVES")
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
                self.shape_functions_derivatives[el.Id-1, count_node, :] = global_gradients_at_nodes[count_node,:]
                self.element_nodes_ids[el.Id-1, count_node] = node.Id-1
                count_node = count_node+1

    def _GetShapeFunctionsDerivatives(self, element):
        gradients = np.asarray(element.GetGeometry().ShapeFunctionDerivatives(1,0))
        jacobian  = np.asarray(element.GetGeometry().Jacobian(0))
        inv_jacobian = np.linalg.inv(jacobian)
        global_gradient = gradients @ inv_jacobian
        return global_gradient 
    
    def _InitializeTopologyOptimizationStepPhysicsSolution(self, first_iteration=False):
        """
        This method Re-Initializes the physical values in order to be consistent with the solution of the current optimization step
        """
        print("\n--|INITIALIZE OPTIMIZATION STEP PHYSICS SOLUTION|")
        self._SetTopologyOptimizationStage(0)
        if (not first_iteration):
                self._ReInitializePhysics()
        self._UpdateDesignParameterAndResistance(self.design_parameter)
    
    def _UpdateDesignParameterAndResistance(self, design_parameter):
        self._UpdateDesignParameter(design_parameter)
        self._UpdateResistance()

    def _UpdateDesignParameter(self, design_parameter):
        print("--|" + self.topology_optimization_stage_str + "| UPDATE DESIGN PARAMETER")
        self.design_parameter_base = design_parameter
        self._ApplyDesignParameterDiffusiveFilter(design_parameter)
        self._ApplyDesignParameterProjectiveFilter(self.design_parameter_filtered)
        self.design_parameter = self.design_parameter_projected
        self._UpdateDesignParameterVariable()

    def _UpdateDesignParameterVariable(self):
        mp = self._GetComputingModelPart()
        for node in mp.Nodes:
                node.SetValue(KratosCFD.DESIGN_PARAMETER, self.design_parameter[node.Id-1])

    def _InitializeResistance(self, min_value = 0.0, max_value = 1000.0, interpolation_slope = 1.0):
        print("--|" + self.topology_optimization_stage_str + "| INITIALIZE RESISTANCE")
        self.alpha_min = min_value
        self.alpha_max = max_value
        self.q = interpolation_slope
        self.resistance = np.zeros(self.n_nodes)
        self.resistance_derivative_wrt_design = np.zeros(self.n_nodes)
        mp = self._GetComputingModelPart()
        for node in mp.Nodes:
            node.SetValue(KratosCFD.RESISTANCE, 0.0)
    
    def _UpdateResistance(self):
        """
        This method will handle the resistance update, that will be based on the value of the design parameter.
        Now it is used to "play" with the resistance
        """
        print("--|" + self.topology_optimization_stage_str + "| UPDATE RESISTANCE")
        self.resistance, self.resistance_derivative_wrt_design_base = self._ComputeResistance(self.design_parameter)
        self._UpdateResistanceDesignDerivative()
        self._UpdateResistanceVariable()

    def _ComputeResistance(self, design_parameter):
        return self._ComputeConvexResistance(design_parameter)

    def _ComputeConvexResistance(self, design_parameter):
        # resistance = a_min + (a_max-a_min)*q*design_parameter/(q+1-design_parameter)
        resistance = self.alpha_min + (self.alpha_max-self.alpha_min)*(self.q*design_parameter)/(self.q+1-design_parameter)
        resistance_derivative_wrt_design_base = (self.alpha_max-self.alpha_min)*(self.q*(self.q+1))/((self.q+1-design_parameter)**2)
        return resistance, resistance_derivative_wrt_design_base
    
    def _UpdateResistanceDesignDerivative(self):
        mask = self._GetOptimizationDomainNodesMask()
        resistance_derivative_wrt_design_projected = self.resistance_derivative_wrt_design_base * self.design_parameter_projected_derivatives
        self.resistance_derivative_wrt_design = resistance_derivative_wrt_design_projected
        self.resistance_derivative_wrt_design[mask] = ((resistance_derivative_wrt_design_projected.T @ self.nodes_connectivity_matrix).T)[mask]
    
    def _UpdateResistanceVariable(self):
        mp = self._GetComputingModelPart()
        for node in mp.Nodes:
                node.SetValue(KratosCFD.RESISTANCE, self.resistance[node.Id-1])

    def _SolveNaviersStokesProblem(self):  
        """
        This method executes a single NS solution of the Topology Optimization problem loop
        """
        self._RunStageSolutionLoop(1) # NAVIER-STOKES PROBLEM SOLUTION  
    
    def _SolveAdjointNaviersStokesProblem(self):  
        """
        This method executes a single ADJ_NS solution of the Topology Optimization problem loop
        """
        self._RunStageSolutionLoop(2) # ADJOINT NAVIER-STOKES PROBLEM SOLUTION

    def _RunStageSolutionLoop(self, problem_stage):
        """
        This method executes a single physics solution of the Topology Optimization problem loop
        N.B.: must be called after the creation of the fluid model part
        problem_stage = 1: Navier Stokes solution
        problem_stage = 2: Adjoint Navier-Stokes solution
        """
        if ((problem_stage == 1)): 
            print("\n--|NAVIER-STOKES SOLUTION|")
        elif (problem_stage == 2):
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
            if (self._GetSolver().IsAdjointNavierStokes()):
                self.OutputSolutionStep()
        print("--|" + top_opt_stage_str + "| END SOLUTION LOOP")

    def KeepAdvancingSolutionLoop(self):
        """This method specifies the stopping criteria for breaking the solution loop
        It can be overridden by derived classes
        """
        if self.IsNavierStokesStage(): # NS
            return self.time < self.end_time
        elif self.IsAdjointNavierStokesStage(): #ADJ
            return self.time > self.start_time
        else:
            print("--|" + self.topology_optimization_stage_str + "| Time check outside NS or ADJ_NS solution")
            return False

    def _AdvanceTime(self):
        """
        This methods hadls the time for the Topology Optimization problems, currently its purpos it's just to make thinks work
        """
        time = super()._AdvanceTime()
        if self.IsAdjointNavierStokesStage():
            time = 0.0
        return time

    def _ReInitializePhysics(self):
        """
        This method Re-Initializes the physical values in order to be consistent with the solution of the current optimization step
        """
        print("--|" + self.topology_optimization_stage_str + "| RE-INITIALIZE PHYSICS")
        self._GetComputingModelPart().ProcessInfo.SetValue(Kratos.TIME, 0.0)

    def _SolveMMA(self, design_parameter, n_opt_variables, n_opt_constraints, min_value = 0.0, max_value = 1.0, max_outer_it = 1, kkt_tolerance = 1e-12):
        print("--|" + self.topology_optimization_stage_str + "| SOLVE MMA")
        # MMA PARAMETERS INITIALIZATION
        n = n_opt_variables
        m = n_opt_constraints
        eeen = np.ones((n, 1))
        eeem = np.ones((m, 1))
        zeron = np.zeros((n, 1))
        zerom = np.zeros((m, 1))
        xval = design_parameter.reshape(-1, 1)
        xold1 = xval.copy()
        xold2 = xval.copy()
        xmin = min_value * eeen
        xmax = max_value * eeen
        low = xmin.copy()
        upp = xmax.copy()
        move = 0.2
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
            f0val, df0dx, fval, dfdx = self._UpdateOptimizationProblem(xval.flatten())
            # outvector1 = np.array([outeriter, innerit, f0val, fval])
            # outvector2 = xval.flatten()
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
            f0val, df0dx, fval, dfdx = self._UpdateOptimizationProblem(xval.flatten())
            # The residual vector of the KKT conditions is calculated
            residu, kktnorm, residumax = MMA.kktcheck(
                m, n, xmma, ymma, zmma, lam, xsi, eta, mu, zet, s, xmin, xmax, df0dx, fval, dfdx, a0, a, c, d)
            # outvector1 = np.array([outeriter, innerit, f0val, fval])
            # outvector2 = xval.flatten()
        new_design_parameter = self._InsertDesignParameterFromOptimizationDomain(xval.flatten())
        self._UpdateDesignParameterAndResistance(new_design_parameter)

    def _EvaluateOptimizationProblem(self, design_parameter = [], first_iteration = False, print_results = False):
        print("\n--|EVALUATE OPTIMIZATION PROBLEM|")
        if (len(design_parameter) != 0):
            self._UpdateDesignParameterAndResistance(design_parameter)
        self._EvaluateFunctionalAndDerivatives(first_iteration, print_results)
        self._EvaluateConstraintsAndDerivatives()
        if (print_results):
            self._PrintOptimizationProblem()
    
    def _UpdateOptimizationProblem(self, optimization_domain_design_parameter):
        design_parameter = self._InsertDesignParameterFromOptimizationDomain(optimization_domain_design_parameter)
        self._EvaluateOptimizationProblem(design_parameter, print_results=False)
        return self.functional, self._ExtractVariableInOptimizationDomain(self.functional_derivatives_wrt_design), self.constraints, self.constraints_derivatives_wrt_design
    
    def _EvaluateConstraintsAndDerivatives(self):
         self._EvaluateVolumeConstraintAndDerivative()

    def _EvaluateFunctionalAndDerivatives(self, first_iteration=False, print_functional=False):
        """
        This method is used to evaluate the functional value and derivatives w.r.t. the design parameter
        """
        self._SetTopologyOptimizationStage(3)
        self._EvaluateFunctional(first_iteration, print_functional)
        self._EvaluateFunctionalDerivatives()
    
    def _EvaluateFunctional(self, first_iteration=False, print_functional=False):
        """
        This method is used to evaluate the functional value
        # Functionals Database
        # 0: resistance  : int_{\Omega}{alpha*||u||^2}
        # 1: strain-rate : int_{\Omega}{2*mu*||S||^2} , with S = 1/2*(grad(u)+grad(u)^T) strain-rate tensor
        # 2: vorticity   : int_{\Omega}{2*mu*||R||^2} = int_{\Omega}{mu*||curl(u)||^2} , curl(u) = vorticity vector, R = 1/2*(grad(u)-grad(u)^T) rotation-rate tensor
        # 3: ...
        # 4: ...
        # 5: ...
        # 5 is just a number big enough to contain the acutal database of functionals
        """
        self._SetTopologyOptimizationStage(3)
        print("--|" + self.topology_optimization_stage_str + "| EVALUATE FUNCTIONAL VALUE")
        mp = self._GetComputingModelPart()
        velocity = np.asarray(Kratos.VariableUtils().GetSolutionStepValuesVector(mp.Nodes, Kratos.VELOCITY, 0, self.dim)).reshape(self.n_nodes, self.dim)
        if (abs(self.functional_weights[0]) > 1e-10):
            self._EvaluateResistanceFunctional(velocity, first_iteration, print_functional)
        if (abs(self.functional_weights[1]) > 1e-10):
            self._EvaluateStrainRateFunctional(velocity, first_iteration, print_functional)
        if (abs(self.functional_weights[2]) > 1e-10):
            self._EvaluateVorticityFunctional(velocity, first_iteration, print_functional)
        self.functional = np.dot(self.functional_weights, self.functionals)
        self.weighted_functionals = self.functional_weights * self.functionals
        if (first_iteration):
            self.initial_functional_abs_value = abs(self.functional)
            if (abs(self.initial_functional_abs_value) < 1e-10):
                self.initial_functional_abs_value = 1.0
            else:
                self.initial_functional_abs_value = 1.0
            self.initial_functionals_abs_value = self.functionals
            self.initial_weighted_functionals_abs_value = self.weighted_functionals
        self.functional = self.functional / self.initial_functional_abs_value

    def _EvaluateResistanceFunctional(self, velocity, first_iteration=False, print_functional=False):
        """
        This method computes the resistance functional: int_{\Omega}{\\alpha||u||^2}
        """
        nodal_velocity_norm = np.linalg.norm(velocity, axis=1)
        integrand = self.resistance * (nodal_velocity_norm**2) #component-wise multiplication
        self.functionals[0] = np.dot(self.nodal_domain_sizes, integrand)
        if (first_iteration):
            self.initial_functionals_values[0] = self.functionals[0] 
        if (print_functional):
            print("--|" + self.topology_optimization_stage_str + "| ---> Resistance Functional (no weight):", self.functionals[0])
        else:
            print("--|" + self.topology_optimization_stage_str + "| --> Resistance Functional")

    def _EvaluateStrainRateFunctional(self, velocity, first_iteration=False, print_functional=False):
        """
        This method computes the Strain-Rate functional: int_{\Omega}{\\2*mu*||1/2*[grad(u)+grad(u)^T]||^2}
        """
        vel = velocity[self.element_nodes_ids[:]]
        vel_gradient = np.matmul(np.transpose(vel, axes=(0,2,1)), self.shape_functions_derivatives)
        vel_symmetric_gradient = 1.0/2.0 * (vel_gradient+(np.transpose(vel_gradient, axes=(0,2,1))))
        vel_symmetric_gradient_norm_squared = (np.linalg.norm(vel_symmetric_gradient, ord='fro', axis=(1, 2)))**2
        mu = self._GetComputingModelPart().Elements[1].Properties.GetValue(Kratos.DYNAMIC_VISCOSITY)
        self.functionals[1] = 2*mu* np.dot(vel_symmetric_gradient_norm_squared, self.elemental_domain_size)
        if (first_iteration):
            self.initial_functionals_values[1] = self.functionals[1]
        if (print_functional):
            print("--|" + self.topology_optimization_stage_str + "| ---> Strain-Rate Functional (no weight):", self.functionals[1])
        else:
            print("--|" + self.topology_optimization_stage_str + "| --> Strain-Rate Functional")

    def _EvaluateVorticityFunctional(self, velocity, first_iteration=False, print_functional=False):
        """
        This method computes the Vorticity functional: int_{\Omega}{\\2*mu*||1/2*[grad(u)-grad(u)^T]||^2}
        """
        vel = velocity[self.element_nodes_ids[:]]
        vel_gradient = np.matmul(np.transpose(vel, axes=(0,2,1)), self.shape_functions_derivatives)
        vel_antisymmetric_gradient = 1.0/2.0 * (vel_gradient-(np.transpose(vel_gradient, axes=(0,2,1))))
        vel_antisymmetric_gradient_norm_squared = (np.linalg.norm(vel_antisymmetric_gradient, ord='fro', axis=(1, 2)))**2
        mu = self._GetComputingModelPart().Elements[1].Properties.GetValue(Kratos.DYNAMIC_VISCOSITY)
        self.functionals[2] = 2*mu* np.dot(vel_antisymmetric_gradient_norm_squared, self.elemental_domain_size)
        if (first_iteration):
            self.initial_functionals_values[2] = self.functionals[2]
        if (print_functional):
            print("--|" + self.topology_optimization_stage_str + "| ---> Vorticity Functional: (no weight)", self.functionals[2])    
        else:
            print("--|" + self.topology_optimization_stage_str + "| --> Vorticity Functional")
    
    def _EvaluateFunctionalDerivatives(self, first_it = False):
        """
        This method is used to evaluate the functional derivatives w.r.t the design parameter
        """
        self._SetTopologyOptimizationStage(3)
        print("--|" + self.topology_optimization_stage_str + "| EVALUATE FUNCTIONAL DERIVATIVES")
        mp = self._GetComputingModelPart()
        velocity = np.asarray(Kratos.VariableUtils().GetSolutionStepValuesVector(mp.Nodes, Kratos.VELOCITY, 0, self.dim)).reshape(self.n_nodes, self.dim)
        velocity_adjoint= np.asarray(Kratos.VariableUtils().GetSolutionStepValuesVector(mp.Nodes, KratosCFD.VELOCITY_ADJ, 0, self.dim)).reshape(self.n_nodes, self.dim)
        # NEXT LINE EXPLANATION:
        # df/dgamma_i = f_weight[0]*dAlpha/dgamma[i]*||vel[i]||^2 + dAlpha/dgamma[i]*(vel[i].dot(vel_adj[i])) = dAlpha/dgamma[i]*vel[i].dot(f_weight[0]*vel[i] + vel_adj[i])
        # then component by component is multiplied by the nodal domain_size of influence
        # !!!!!!!! ASK RICCARDO IF IT IS CORRECT !!!!!!!!!
        self.functional_derivatives_wrt_design = self.resistance_derivative_wrt_design * np.sum(velocity * (self.functional_weights[0]*velocity + velocity_adjoint), axis=1) * self.nodal_domain_sizes 
        self.functional_derivatives_wrt_design = self.functional_derivatives_wrt_design / self.initial_functional_abs_value
        for node in mp.Nodes:
            node.SetValue(KratosCFD.FUNCTIONAL_DERIVATIVE, self.functional_derivatives_wrt_design[node.Id-1])
        if (first_it):
            self.initial_functional_derivatives_wrt_design = self.functional_derivatives_wrt_design
        self.functional_derivatives_wrt_design = np.asarray([self.functional_derivatives_wrt_design]).T

    def _EvaluateVolumeConstraintAndDerivative(self, ):
        self.volume_fraction = 1.0 - np.dot(self.design_parameter, self.nodal_optimization_domain_sizes)/self.optimization_domain_size
        self.volume_constraint = self.volume_fraction - self.max_volume_fraction
        nodal_optimization_domain_sizes_extracted = self._ExtractVariableInOptimizationDomain(self.nodal_optimization_domain_sizes)
        self.volume_constraint_derivatives_wrt_design = -1.0 * nodal_optimization_domain_sizes_extracted / self.optimization_domain_size
        self.constraints[0] = self.volume_constraint
        self.constraints_derivatives_wrt_design[0,:] = self.volume_constraint_derivatives_wrt_design

    def _EvaluateWSSConstraintAndDerivative(self):
        dumb = 0
        return dumb
    
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
            self.topology_optimization_stage_str = "N-S"
        elif self.CheckOptimizationStage(2):
            self.topology_optimization_stage_str = "ADJ"
        elif self.CheckOptimizationStage(3):
            self.topology_optimization_stage_str = "OPT"
        else:
            self.topology_optimization_stage_str = "ERROR"
        self._GetComputingModelPart().ProcessInfo.SetValue(KratosCFD.FLUID_TOP_OPT_PROBLEM_STAGE, self.topology_optimization_stage)

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
    
    def IsInitializeStage(self):
        """
        This method returns TRUE the topology optimization stage is 0 <=> INIT
        """ 
        return self.CheckOptimizationStage(0)
    
    def IsNavierStokesStage(self):
        """
        This method returns TRUE the topology optimization stage is 1 <=> NS
        """
        return self.CheckOptimizationStage(1)

    def IsAdjointNavierStokesStage(self):
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
    
    def _GetSubModelPart(self, sub_model_part_name):
        """
        This method returns the selected sub model part for the current solver
        """  
        if (self._GetMainModelPart().HasSubModelPart(sub_model_part_name)):  
            return self._GetMainModelPart().GetSubModelPart(sub_model_part_name)
        else:
            return None
    
    def _GetOptimizationDomain(self):
        """
        This method returns the optimization domain model part for the current solver
        """  
        return self._GetSubModelPart("GENERIC_domain-optimization_domain")
    
    def _GetNonOptimizationDomain(self):
        """
        This method returns the non-optimization domain model part for the current solver
        """  
        return self._GetSubModelPart("GENERIC_domain-non_optimization_domain")

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
        print("\n--|PRINT OPTIMIZATION PROBLEM DATA|")
        self._PrintFunctionals()
        self._PrintConstraints()
    
    def _PrintFunctionals(self):
        print("--|" + self.topology_optimization_stage_str + "| TOTAL FUNCTIONAL  :", self.functional)
        print("--|" + self.topology_optimization_stage_str + "| INITIAL FUNCTIONAL:", self.initial_functional_abs_value)
        if (abs(self.functional_weights[0]) > 1e-10):
            print("--|" + self.topology_optimization_stage_str + "| ---> Resistance Functional (" + str(self.functional_weights[0]) + "):", self.weighted_functionals[0]/self.initial_functional_abs_value)
        if (abs(self.functional_weights[1]) > 1e-10):
            print("--|" + self.topology_optimization_stage_str + "| ---> Strain-Rate Functional (" + str(self.functional_weights[1]) + "):", self.weighted_functionals[1]/self.initial_functional_abs_value)
        if (abs(self.functional_weights[2]) > 1e-10):
            print("--|" + self.topology_optimization_stage_str + "| ---> Vorticity Functional (" + str(self.functional_weights[2]) + "):", self.weighted_functionals[2]/self.initial_functional_abs_value)
        
    def _PrintFunctionalsToFile(self):
        with open("functional_history.txt", "a") as file:
            file.write(str(self.functional) + " ")
            for ifunc in range(self.n_functionals): 
                file.write(str(self.weighted_functionals[ifunc]) + " ")
            file.write("\n")

    def _EvaluateDesignParameterChange(self):
        old_design_parameter_norm = np.linalg.norm(1-self.old_design_parameter)
        if (old_design_parameter_norm > 1e-10):
            # change evaluated on the amount of fluid, that's why we have ||(1-design)-(1-old_design)|| / ||1-old_design||
            self.design_parameter_change = np.linalg.norm(self.design_parameter-self.old_design_parameter)
            self.design_parameter_change /= old_design_parameter_norm
        else:
            pass
        if (self.opt_it > 1):
            design_parameter_converged = (self.design_parameter_change < self.design_parameter_change_toll)
        else:
            design_parameter_converged = False
        print("--|" + self.topology_optimization_stage_str + "| DESIGN PARAMETER CHANGE:", self.design_parameter_change)
        return design_parameter_converged        

    def _ResetFunctionalOutput(self):
        with open("functional_history.txt", "w") as file:
            file.write("TOTAL | RESISTANCE | STRAIN-RATE | VORTICITY | NONE | NONE\n")
            file.write("1 ")
            for ifunc in range(self.n_functionals):
                file.write(str(self.functional_weights[ifunc]) + " ")
            file.write("\n")

    def _PrintConstraints(self):
        self._PrintVolumeConstraint()
    
    def _PrintVolumeConstraint(self):
        print("--|" + self.topology_optimization_stage_str + "| VOLUME FRACTION:", self.volume_fraction)
        print("--|" + self.topology_optimization_stage_str + "| ---> Volume Constraint:", self.volume_constraint)

    def PrintAnalysisStageProgressInformation(self):
        Kratos.Logger.PrintInfo(self._GetSimulationName(), "TOTAL STEP: ", self._GetComputingModelPart().ProcessInfo[Kratos.STEP])
        if self.IsNavierStokesStage(): # NS
            Kratos.Logger.PrintInfo(self._GetSimulationName(), "NS STEP: ", self._GetComputingModelPart().ProcessInfo[KratosCFD.FLUID_TOP_OPT_NS_STEP])
        elif self.IsAdjointNavierStokesStage(): # ADJ
            Kratos.Logger.PrintInfo(self._GetSimulationName(), "ADJ NS STEP: ", self._GetComputingModelPart().ProcessInfo[KratosCFD.FLUID_TOP_OPT_ADJ_NS_STEP])
        else:
            Kratos.Logger.PrintInfo(self._GetSimulationName(), "Invalid Value of the Topology Optimization Stage. TOP_OPT_STAGE: ", self.topology_optimization_stage)
        Kratos.Logger.PrintInfo(self._GetSimulationName(), "TIME: ", self.time)

    def _ComputeDesignParameterDiffusiveFilterUtilities(self):
        self._BuildNodesConnectivity()
        self.design_parameter_filtered = np.zeros(self.n_nodes)

    def _ComputeDesignParameterProjectiveFilterUtilities(self):
        self.design_parameter_projected = np.zeros(self.n_nodes)
        self.design_parameter_projected_derivatives = np.ones(self.n_nodes)

    # def _BuildNodesConnectivity(self, r_max = 0.0):
    #     self.diffusive_filter_radius = r_max
    #     self.nodes_connectivity = [[i] for i in range(self.n_nodes)]
    #     self.nodes_connectivity_weights = [[self.diffusive_filter_radius] for i in range(self.n_nodes)]
    #     self.nodes_connectivity_weights_as_nb = [[self.diffusive_filter_radius] for i in range(self.n_nodes)]
    #     self.nodes_connectivity_weights_sum = np.zeros(self.n_nodes)
    #     mp = self._GetComputingModelPart()
    #     for elem in mp.Elements:
    #         elem_geometry = elem.GetGeometry()
    #         for node_i in elem_geometry:
    #             node_i_coords = self._GetNodeCoordinates(node_i)
    #             for node_j in elem_geometry:
    #                 if (node_j.Id != node_i.Id):
    #                     node_j_coords = self._GetNodeCoordinates(node_j)
    #                     self.nodes_connectivity[node_i.Id-1].append(node_j.Id-1)
    #                     distance = np.linalg.norm(node_i_coords-node_j_coords)
    #                     nodal_weight = max(0.0, self.diffusive_filter_radius-distance)
    #                     self.nodes_connectivity_weights[node_i.Id-1].append(nodal_weight)
    #                     self.nodes_connectivity_weights_as_nb[node_i.Id-1].append(nodal_weight)

    #     # normalize nodal weights
    #     for node in mp.Nodes:
    #         self.nodes_connectivity_weights[node.Id-1] = np.asarray(self.nodes_connectivity_weights[node.Id-1])
    #         nodal_weigths_sum = np.sum(self.nodes_connectivity_weights[node.Id-1])
    #         self.nodes_connectivity_weights_sum[node.Id-1] = nodal_weigths_sum
    #         self.nodes_connectivity_weights[node.Id-1] /= nodal_weigths_sum

    #     # normalize nodeal weights ad neighbours (the weight as nb is not obtained normalizing by the ccurrent node weiights sum, but eveery weight must be divivded by the sum of the relative nodes weights)
    #     for inode in range(self.n_nodes):
    #         self.nodes_connectivity_weights_as_nb[inode] = np.asarray(self.nodes_connectivity_weights_as_nb[inode])
    #         for jnode in range(len(self.nodes_connectivity[inode])):
    #             jnode_id = self.nodes_connectivity[inode][jnode]
    #             jnodal_weights_sum = self.nodes_connectivity_weights_sum[jnode_id]
    #             self.nodes_connectivity_weights_as_nb[inode][jnode] /= jnodal_weights_sum
                
    # def _ApplyDesignParameterDiffusiveFilter(self, design_parameter):
    #     opt_mp = self._GetOptimizationDomain()
    #     for node in opt_mp.Nodes:
    #         if (self.is_only_opt_domain_node(node.Id-1)):
    #             temp_nb = self.nodes_connectivity[node.Id-1]
    #             temp_weights = self.nodes_connectivity_weights[node.Id-1]
    #             temp_design_parameter = design_parameter[temp_nb]
    #             self.design_parameter_filtered[node.Id-1] = np.dot(temp_weights, temp_design_parameter)

    def _BuildNodesConnectivity(self, radius = 1e-10):
        print("--|" + self.topology_optimization_stage_str + "| ---> Build Nodes Connectivity Matrix")
        mp = self._GetComputingModelPart()
        points = self._GetNodesSetCoordinates(mp.Nodes)
        nodes_tree = KDTree(points)
        self.nodes_connectivity_matrix = nodes_tree.sparse_distance_matrix(nodes_tree, radius, output_type="dok_matrix").tocsr()
        self.nodes_connectivity_matrix *= -1
        self.nodes_connectivity_matrix.data += radius
        self.nodes_connectivity_weigths_sum = np.array(self.nodes_connectivity_matrix.sum(axis=1)).flatten()  # Sum of each row as a 1D array
        
        # Normalization step: Divide non-zero entries by the corresponding row sum
        for i in range(self.nodes_connectivity_matrix.shape[0]):  # Iterate over rows
            start = self.nodes_connectivity_matrix.indptr[i]
            end = self.nodes_connectivity_matrix.indptr[i + 1]
            self.nodes_connectivity_matrix.data[start:end] /= self.nodes_connectivity_weigths_sum[i]  # Normalize non-zero entries of the row

    def _ApplyDesignParameterDiffusiveFilter(self, design_parameter):
        print("--|" + self.topology_optimization_stage_str + "| --> Apply Diffusive Filter")
        mask = self._GetOptimizationDomainNodesMask()
        self.design_parameter_filtered = design_parameter
        self.design_parameter_filtered[mask] = (self.nodes_connectivity_matrix @ design_parameter)[mask]

    def _ApplyDesignParameterProjectiveFilter(self, design_parameter, mean_min = 0.2, mean_max = 0.5, projection_slope = 2):
        # DO NOTHING
        print("--|" + self.topology_optimization_stage_str + "| --> Apply Projective Filter")
        mask = self._GetOptimizationDomainNodesMask()
        self.design_parameter_projected = design_parameter
        design_parameter_mean = min(max(mean_min,1.0-self.volume_fraction),mean_max)
        factor = np.tanh(projection_slope*design_parameter_mean) + np.tanh(projection_slope*(1-design_parameter_mean))
        constant = np.tanh(projection_slope*design_parameter_mean)
        projection_argument = projection_slope*(design_parameter-design_parameter_mean)
        self.design_parameter_projected[mask] = ((constant + np.tanh(projection_argument)) / factor)[mask]
        self.design_parameter_projected_derivatives[mask] = ((projection_slope / factor) / (np.cosh(projection_argument)**2))[mask]
                


    













        