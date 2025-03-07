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
# Import Kratos Applications
import KratosMultiphysics.FluidDynamicsApplication as KratosCFD
import KratosMultiphysics.MeshingApplication as KratosMMG
# Import Kratos Analysis and Solvers
from KratosMultiphysics.FluidDynamicsApplication.fluid_dynamics_analysis import FluidDynamicsAnalysis
from KratosMultiphysics.FluidDynamicsApplication import fluid_topology_optimization_solver
# Import Kratos Processes
from KratosMultiphysics import ComputeNodalGradientProcess

class FluidTopologyOptimizationAnalysis(FluidDynamicsAnalysis):
    def __init__(self,model,parameters):
        self.topology_optimization_stage = 0
        self.topology_optimization_stage_str = "INIT"
        super().__init__(model,parameters) 
        # self._CreateTopologyOptimizationSolvers() # currently it is a useless method 
        self._SetMinMaxIt()  
        self._SetTopologyOptimizationName()

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
        # # Set the Functionals Weigths for the Optimization
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

    def _SetFunctionalWeights(self, weights = [1, 1, 0, 0, 0, 0, 0, 0, 0]):
        # set future transport functinals to zero
        weights[3]=0.0
        weights[4]=0.0
        weights[5]=0.0
        weights[6]=0.0
        weights[7]=0.0
        weights[8]=0.0
        self.n_functionals = len(weights)
        self.initial_functionals_values = np.zeros(self.n_functionals)
        self.functionals = np.zeros(self.n_functionals)
        self.functional_weights = self._NormalizeFunctionalWeights(np.asarray(weights))
        self._GetComputingModelPart().ProcessInfo.SetValue(KratosMultiphysics.FUNCTIONAL_WEIGHTS, self.functional_weights)

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
        print("\n------------------------------------------------------------------------------")
        print(  "--|", self.topology_optimization_name, "TOPOLOGY OPTIMIZATION PREPROCESSING")
        print("------------------------------------------------------------------------------")
        print("--|INITIALIZE|")
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
            print("\n------------------------------------------------------------------------------")
            print(  "--|", self.topology_optimization_name, "TOPOLOGY OPTIMIZATION SOLUTION LOOP. IT:", self.opt_it)
            print("------------------------------------------------------------------------------")
            self.old_design_parameter = self.design_parameter
            self._SolveOptimizer()
            self._SolveTopologyOptimizationStepPhysics()
            self._EvaluateOptimizationProblem(print_results=True)
            end_solution = self._IsTopologyOptimizationSolutionEnd()
            if (end_solution):
                print("\n------------------------------------------------------------------------------")
                print("--| ENDING", self.topology_optimization_name, "TOPOLOGY OPTIMIZATION SOLUTION LOOP")
                if (self.converged):
                    print("--| ---> CONVERGED!")
                else:
                    print("--| ---> Reached Max Number of Iterations")
                self._Remesh()
                print("------------------------------------------------------------------------------\n")
            self._PrintSolution()
            self.first_iteration = False
            
    def _IsTopologyOptimizationSolutionEnd(self):
        design_parameter_converged = self._EvaluateDesignParameterChange()
        volume_constraint_valid = not (self.volume_constraint > 0.0)
        self.converged = design_parameter_converged and volume_constraint_valid
        return not (((self.opt_it < self.max_it) and (not self.converged)) or (self.opt_it < self.min_it))

    def _SolveOptimizer(self):
        print("\n--|OPTIMIZER|")
        self._SetTopologyOptimizationStage(3)
        if (self.first_iteration):
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
        self._InitializeTopologyOptimizationStepPhysicsSolution()
        self._CheckMaterialProperties()
        self._UpdateRelevantPhysicsVariables()
        self._SolvePhysicsProblem() # NAVIER-STOKES PROBLEM SOLUTION
        if (self.first_iteration):
            self._SetFunctionalWeights()
            self._ResetFunctionalOutput()
        self._UpdateRelevantAdjointVariables()
        self._SolveAdjointProblem() # ADJOINT NAVIER-STOKES PROBLEM SOLUTION   
    
    def _GeometricalPreprocessing(self):
        """
        This method preprocess all the useful values and quantities for the topology optimization solution
        """
        print("--|" + self.topology_optimization_stage_str + "| GEOMETRICAL PREPROCESSING")
        mp = self._GetComputingModelPart()
        self.dim = mp.ProcessInfo.GetValue(KratosMultiphysics.DOMAIN_SIZE)
        self.nodes_in_element = mp.Elements[1].GetGeometry().PointsNumber()
        self._OrderNodes()
        self._OrderElements()
        self._ComputeNodalDomainSizes()
    
    def _InitializeOptimization(self):  
        print("--|" + self.topology_optimization_stage_str + "| INITIALIZE OPTIMIZATION")
        self._OptimizationGeometricalPreprocessing()
        self.opt_it = 0
        self._SetDesignParameterChangeTolerance()
        self.n_optimization_constraints = 1  
        self._InitializeConstraints()
        self._InitializeRemeshing()
    
    def _SetDesignParameterChangeTolerance(self, design_parameter_change_toll = 0.0):
        self.design_parameter_change_toll = design_parameter_change_toll

    def _ResetConstraints(self):
        self.constraints = np.zeros(self.n_optimization_constraints)  
        self.constraints_derivatives_wrt_design = np.zeros((self.n_optimization_constraints, self.n_opt_design_parameters))

    def _InitializeConstraints(self):
        self._ResetConstraints()
        self._InitializeVolumeConstraint()

    def _InitializeVolumeConstraint(self, is_fluid_volume=True):
        self.is_fluid_volume_constraint = is_fluid_volume

    def _InitializeDomainDesign(self):
        print("--|" + self.topology_optimization_stage_str + "| INITIALIZE DOMAIN DESIGN")
        self._SetMaxDomainVolumeFraction()
        self._InitializeDomainDesignParameter()
        self._InitializePhysicsParameters()

    def _InitializePhysicsParameters(self, resistance_parameters=[1e4, 0.0, 1.0]):
        self._SetPhysicsParametersInterpolationMethods()
        self._InitializeResistance(resistance_parameters)

    def _SetPhysicsParametersInterpolationMethods(self, resistance="hyperbolic"):
        self._SetResistanceInterpolationMethod(resistance)

    def _SetResistanceInterpolationMethod(self, interpolation_method):
        self.resistance_interpolation_method = interpolation_method
    
    def _SetMaxDomainVolumeFraction(self, max_volume_fraction = 1.0):
        self.max_volume_fraction = max_volume_fraction
        # self.volume_fraction = 0.0

    def _InitializeDomainDesignParameter(self, initial_value = 0.0):
        """
        This method will handle the design parameter initialization across the whole domain.
        self.design_parameter is defined in every node of the mesh.
        The optimization process will update only the value of the nodes belonging to the optimization_domain sub model part
        """
        mask = self._GetOptimizationDomainNodesMask()
        self.volume_fraction = initial_value
        self.design_parameter = np.zeros(self.n_nodes) 
        self.design_parameter[mask] = initial_value
        mp = self._GetComputingModelPart()
        for node in mp.Nodes:
            design = self.design_parameter[node.Id-1]
            node.SetValue(KratosMultiphysics.DESIGN_PARAMETER, design)
            distance = design-self.remeshing_levelset
            node.SetSolutionStepValue(KratosMultiphysics.DISTANCE, distance)

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
        This method compute the nodal domain size - not vectorized but it is done once
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
        self._UpdateNodalAreaVariable()
        
    def _UpdateNodalAreaVariable(self):
        mp = self._GetComputingModelPart()
        for node in mp.Nodes:
            node.SetValue(KratosMultiphysics.NODAL_AREA, self.nodal_domain_sizes[node.Id-1])

    def _ComputeScalarVariableNodalGradient(self, scalar_variable, gradient_variable):
        mp = self._GetComputingModelPart()
        gradient_process = ComputeNodalGradientProcess(mp, scalar_variable, gradient_variable, KratosMultiphysics.NODAL_AREA)
        gradient_process.Execute()

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
    
    def _InitializeTopologyOptimizationStepPhysicsSolution(self):
        """
        This method Re-Initializes the physical values in order to be consistent with the solution of the current optimization step
        """
        print("\n--|INITIALIZE OPTIMIZATION STEP PHYSICS SOLUTION|")
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
        print("--|" + self.topology_optimization_stage_str + "| UPDATE DESIGN PARAMETER")
        self.design_parameter_base = design_parameter
        self._ApplyDesignParameterDiffusiveFilter(design_parameter)
        self._ApplyDesignParameterProjectiveFilter(self.design_parameter_filtered)
        self.design_parameter = self.design_parameter_projected
        self._UpdateDesignParameterVariable()

    def _UpdateDesignParameterVariable(self):
        mp = self._GetComputingModelPart()
        for node in mp.Nodes:
            design = self.design_parameter[node.Id-1]
            node.SetValue(KratosMultiphysics.DESIGN_PARAMETER, design)
            distance = design-self.remeshing_levelset
            node.SetSolutionStepValue(KratosMultiphysics.DISTANCE, distance)

    def _ResetPhysicsParameters(self):
        self._ResetResistance()

    def _InitializeResistance(self, resistance_parameters_values):
        print("--|" + self.topology_optimization_stage_str + "| INITIALIZE RESISTANCE")
        self.resistance_parameters = resistance_parameters_values
        self._ResetResistance()
        self._UpdateResistanceVariable()

    def _ResetResistance(self):
        self.resistance = np.zeros(self.n_nodes)
        self.resistance_derivative_wrt_design = np.zeros(self.n_nodes)
    
    def _UpdateResistance(self):
        """
        This method handles the resistance update.
        """
        print("--|" + self.topology_optimization_stage_str + "| UPDATE RESISTANCE")
        self.resistance, self.resistance_derivative_wrt_design_base = self._ComputeResistance(self.design_parameter)
        self._UpdateResistanceDesignDerivative()
        self._UpdateResistanceVariable()

    def _ComputeResistance(self, design_parameter):
        return self._ComputePhysicsParameter(self.resistance_interpolation_method, self.resistance_parameters, design_parameter)
    
    def _ComputePhysicsParameter(self, interpolation_method, physics_parameter_values, design_parameter):
        if (interpolation_method == ("POLYNOMIAL").lower()):
            return self._ComputePolynomialPhysicsParameter(physics_parameter_values, design_parameter)
        else:
            return self._ComputeHyperbolicPhysicsParameter(physics_parameter_values, design_parameter)
        
    def _ComputeHyperbolicPhysicsParameter(self, physics_parameter_values, design_parameter):
        value_void = physics_parameter_values[0]
        value_solid = physics_parameter_values[1]
        slope       = physics_parameter_values[2]
        if (value_void <= value_solid): 
            physics_parameter = value_void + (value_solid-value_void)*(slope*design_parameter)/(slope+1-design_parameter)
            physics_parameter_derivative_wrt_design_base = (value_solid-value_void)*(slope*(slope+1))/((slope+1-design_parameter)**2)
        else:
            physics_parameter = value_solid - (value_solid-value_void)*(slope*(1-design_parameter))/(slope+design_parameter)
            physics_parameter_derivative_wrt_design_base = (value_solid-value_void)*(slope*(slope+1))/((slope+design_parameter)**2)
        return physics_parameter, physics_parameter_derivative_wrt_design_base
    
    def _ComputePolynomialPhysicsParameter(self, physics_parameter_values, design_parameter):
        value_void  = physics_parameter_values[0]
        value_solid = physics_parameter_values[1]
        power       = physics_parameter_values[2]
        if (power < 1.0):
            power = 1.0
        power_multiplier = 10.0
        power_max   = power_multiplier*power
        min_it = 5.0
        max_it = 100.0
        if (self.opt_it < min_it):
            eff_power = power
        elif (self.opt_it <= max_it):
            eff_power = power + (power_max-power)*(self.opt_it-min_it)/(max_it-min_it)
        else:
            eff_power = power_max
        if (value_void <= value_solid): 
            physics_parameter = value_void + (value_solid-value_void)*(design_parameter**eff_power)
            physics_parameter_derivative_wrt_design_base = eff_power*(value_solid-value_void)*(design_parameter**(eff_power-1))
        else:
            physics_parameter = value_solid - (value_solid-value_void)*((1.0-design_parameter)**eff_power)
            physics_parameter_derivative_wrt_design_base = eff_power*(value_solid-value_void)*((1.0-design_parameter)**(eff_power-1))
        return physics_parameter, physics_parameter_derivative_wrt_design_base
    
    def _UpdateResistanceDesignDerivative(self):
        mask = self._GetOptimizationDomainNodesMask()
        resistance_derivative_wrt_design_projected = self.resistance_derivative_wrt_design_base * self.design_parameter_projected_derivatives
        self.resistance_derivative_wrt_design = resistance_derivative_wrt_design_projected
        self.resistance_derivative_wrt_design[mask] = self._ApplyDiffusiveFilterDerivative(resistance_derivative_wrt_design_projected)[mask]
    
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

    def KeepAdvancingSolutionLoop(self):
        """This method specifies the stopping criteria for breaking the solution loop
        It can be overridden by derived classes
        """
        if self.IsPhysicsStage(): # NS
            return self.time < self.end_time
        elif self.IsAdjointStage(): #ADJ
            return self.time > self.start_time
        else:
            print("--|" + self.topology_optimization_stage_str + "| Time check outside NS or ADJ_NS solution")
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
        print("--|" + self.topology_optimization_stage_str + "| RE-INITIALIZE PHYSICS")
        self._GetComputingModelPart().ProcessInfo.SetValue(KratosMultiphysics.TIME, 0.0)

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
        self._UpdateDesignParameterAndPhysicsParameters(new_design_parameter)

    def _EvaluateOptimizationProblem(self, design_parameter = [], print_results = False):
        print("\n--|EVALUATE OPTIMIZATION PROBLEM|")
        if (len(design_parameter) != 0):
            self._UpdateDesignParameterAndPhysicsParameters(design_parameter)
        self._EvaluateFunctionalAndDerivatives(print_results)
        self._EvaluateConstraintsAndDerivatives()
        if (print_results):
            self._PrintOptimizationProblem()
    
    def _UpdateOptimizationProblem(self, optimization_domain_design_parameter):
        design_parameter = self._InsertDesignParameterFromOptimizationDomain(optimization_domain_design_parameter)
        self._EvaluateOptimizationProblem(design_parameter, print_results=False)
        return self.functional, self._ExtractVariableInOptimizationDomain(self.functional_derivatives_wrt_design), self.constraints, self.constraints_derivatives_wrt_design
    
    def _EvaluateConstraintsAndDerivatives(self):
         self._EvaluateVolumeConstraintAndDerivative()

    def _EvaluateFunctionalAndDerivatives(self, print_functional=False):
        """
        This method is used to evaluate the functional value and derivatives w.r.t. the design parameter
        """
        self._SetTopologyOptimizationStage(3)
        self._EvaluateFunctional(print_functional)
        self._EvaluateFunctionalDerivatives()
    
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
        print("--|" + self.topology_optimization_stage_str + "| EVALUATE FUNCTIONAL VALUE")
        self.EvaluateFunctionals(print_functional)
        self.functional = np.dot(self.functional_weights, self.functionals)
        self.weighted_functionals = self.functional_weights * self.functionals
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

    def _EvaluateResistanceFunctional(self, print_functional=False):
        """
        This method computes the resistance functional: int_{\Omega}{\\alpha||u||^2}
        """
        mp = self._GetComputingModelPart()
        velocity = np.asarray(KratosMultiphysics.VariableUtils().GetSolutionStepValuesVector(mp.Nodes, KratosMultiphysics.VELOCITY, 0, self.dim)).reshape(self.n_nodes, self.dim)
        nodal_velocity_norm = np.linalg.norm(velocity, axis=1)
        integrand = self.resistance * (nodal_velocity_norm**2) #component-wise multiplication
        self.functionals[0] = np.dot(self.nodal_domain_sizes, integrand)
        if (self.first_iteration):
            self.initial_functionals_values[0] = self.functionals[0] 
        if (print_functional):
            print("--|" + self.topology_optimization_stage_str + "| ---> Resistance Functional (no weight):", self.functionals[0])
        else:
            print("--|" + self.topology_optimization_stage_str + "| ---> Resistance Functional")

    def _EvaluateStrainRateFunctional(self, print_functional=False):
        """
        This method computes the Strain-Rate functional: int_{\Omega}{\\2*mu*||1/2*[grad(u)+grad(u)^T]||^2}
        """
        mp = self._GetComputingModelPart()
        velocity = np.asarray(KratosMultiphysics.VariableUtils().GetSolutionStepValuesVector(mp.Nodes, KratosMultiphysics.VELOCITY, 0, self.dim)).reshape(self.n_nodes, self.dim)
        vel = velocity[self.element_nodes_ids[:]]
        vel_gradient = np.matmul(np.transpose(vel, axes=(0,2,1)), self.shape_functions_derivatives)
        vel_symmetric_gradient = 1.0/2.0 * (vel_gradient+(np.transpose(vel_gradient, axes=(0,2,1))))
        vel_symmetric_gradient_norm_squared = (np.linalg.norm(vel_symmetric_gradient, ord='fro', axis=(1, 2)))**2
        mu = self._GetComputingModelPart().Elements[1].Properties.GetValue(KratosMultiphysics.DYNAMIC_VISCOSITY)
        self.functionals[1] = 2*mu* np.dot(vel_symmetric_gradient_norm_squared, self.elemental_domain_size)
        if (self.first_iteration):
            self.initial_functionals_values[1] = self.functionals[1]
        if (print_functional):
            print("--|" + self.topology_optimization_stage_str + "| ---> Strain-Rate Functional (no weight):", self.functionals[1])
        else:
            print("--|" + self.topology_optimization_stage_str + "| ---> Strain-Rate Functional")

    def _EvaluateVorticityFunctional(self, print_functional=False):
        """
        This method computes the Vorticity functional: int_{\Omega}{\\2*mu*||1/2*[grad(u)-grad(u)^T]||^2}
        """
        mp = self._GetComputingModelPart()
        velocity = np.asarray(KratosMultiphysics.VariableUtils().GetSolutionStepValuesVector(mp.Nodes, KratosMultiphysics.VELOCITY, 0, self.dim)).reshape(self.n_nodes, self.dim)
        vel = velocity[self.element_nodes_ids[:]]
        vel_gradient = np.matmul(np.transpose(vel, axes=(0,2,1)), self.shape_functions_derivatives)
        vel_antisymmetric_gradient = 1.0/2.0 * (vel_gradient-(np.transpose(vel_gradient, axes=(0,2,1))))
        vel_antisymmetric_gradient_norm_squared = (np.linalg.norm(vel_antisymmetric_gradient, ord='fro', axis=(1, 2)))**2
        mu = self._GetComputingModelPart().Elements[1].Properties.GetValue(KratosMultiphysics.DYNAMIC_VISCOSITY)
        self.functionals[2] = 2*mu* np.dot(vel_antisymmetric_gradient_norm_squared, self.elemental_domain_size)
        if (self.first_iteration):
            self.initial_functionals_values[2] = self.functionals[2]
        if (print_functional):
            print("--|" + self.topology_optimization_stage_str + "| ---> Vorticity Functional: (no weight)", self.functionals[2])    
        else:
            print("--|" + self.topology_optimization_stage_str + "| ---> Vorticity Functional")
    
    def _EvaluateFunctionalDerivatives(self):
        """
        This method is used to evaluate the functional derivatives w.r.t the design parameter
        """
        self._SetTopologyOptimizationStage(3)
        print("--|" + self.topology_optimization_stage_str + "| EVALUATE FUNCTIONAL DERIVATIVES")
        self.functional_derivatives_wrt_design = np.asarray([self._ComputeFunctionalDerivatives()]).T
        self._UpdateFunctionalDerivativesVariable()
        if (self.first_iteration):
            self.initial_functional_derivatives_wrt_design = self.functional_derivatives_wrt_design

    def _ComputeFunctionalDerivatives(self):
        temp_functional_derivatives_wrt_design = np.asarray(self.n_nodes)
        temp_functional_derivatives_wrt_design = self._ComputeFunctionalDerivativesFunctionalContribution()
        temp_functional_derivatives_wrt_design += self._ComputeFunctionalDerivativesPhysicsContribution()
        return temp_functional_derivatives_wrt_design
        
    def _UpdateFunctionalDerivativesVariable(self):
        mp = self._GetComputingModelPart()
        for node in mp.Nodes:
            node.SetValue(KratosMultiphysics.FUNCTIONAL_DERIVATIVE, self.functional_derivatives_wrt_design[node.Id-1][0])

    def _ComputeFunctionalDerivativesFunctionalContribution(self):
        return self._ComputeFunctionalDerivativesFluidFunctionalContribution()
    
    def _ComputeFunctionalDerivativesPhysicsContribution(self):
        return self._ComputeFunctionalDerivativesFluidPhysicsContribution()

    def _ComputeFunctionalDerivativesFluidFunctionalContribution(self):
        mp = self._GetComputingModelPart()
        velocity = np.asarray(KratosMultiphysics.VariableUtils().GetSolutionStepValuesVector(mp.Nodes, KratosMultiphysics.VELOCITY, 0, self.dim)).reshape(self.n_nodes, self.dim)
        return self.functional_weights[0]*self.resistance_derivative_wrt_design * np.sum(velocity*velocity, axis=1) * self.nodal_domain_sizes 

    def _ComputeFunctionalDerivativesFluidPhysicsContribution(self):
        mp = self._GetComputingModelPart()
        velocity = np.asarray(KratosMultiphysics.VariableUtils().GetSolutionStepValuesVector(mp.Nodes, KratosMultiphysics.VELOCITY, 0, self.dim)).reshape(self.n_nodes, self.dim)
        velocity_adjoint= np.asarray(KratosMultiphysics.VariableUtils().GetSolutionStepValuesVector(mp.Nodes, KratosMultiphysics.VELOCITY_ADJ, 0, self.dim)).reshape(self.n_nodes, self.dim)        
        return self.resistance_derivative_wrt_design * np.sum(velocity*velocity_adjoint, axis=1) * self.nodal_domain_sizes

    def _EvaluateVolumeConstraintAndDerivative(self):
        mask = self._GetOptimizationDomainNodesMask()
        if (self.is_fluid_volume_constraint):
            self.volume_fraction = 1.0 - np.dot(self.design_parameter, self.nodal_optimization_domain_sizes)/self.optimization_domain_size
            volume_constraint_derivatives_wrt_design_base = -1.0 * self.nodal_optimization_domain_sizes / self.optimization_domain_size
        else:
            self.volume_fraction = np.dot(self.design_parameter, self.nodal_optimization_domain_sizes)/self.optimization_domain_size
            volume_constraint_derivatives_wrt_design_base = self.nodal_optimization_domain_sizes / self.optimization_domain_size
        self.volume_constraint = self.volume_fraction - self.max_volume_fraction
        volume_constraint_derivatives_wrt_design_projected = volume_constraint_derivatives_wrt_design_base * self.design_parameter_projected_derivatives
        self.constraints[0] = self.volume_constraint
        self.constraints_derivatives_wrt_design[0,:] = self._ApplyDiffusiveFilterDerivative(volume_constraint_derivatives_wrt_design_projected)[mask]

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
        opt_domain = self._GetSubModelPart("GENERIC_domain-optimization_domain")
        if (opt_domain is None):
            return self._GetMainModelPart()
        else:
            return opt_domain
    
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
        print("--|" + self.topology_optimization_stage_str + "| INITIAL FUNCTIONAL:", self.initial_functional)
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
        print("--|" + self.topology_optimization_stage_str + "| DESIGN PARAMETER CHANGE:", self.design_parameter_change)
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
    
    def _PrintVolumeConstraint(self):
        print("--|" + self.topology_optimization_stage_str + "| VOLUME FRACTION:", self.volume_fraction)
        print("--|" + self.topology_optimization_stage_str + "| ---> Volume Constraint:", self.volume_constraint)

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

    def _InitializeDiffusiveFilter(self, apply_diffusive_filter=True, diffusive_filter_radius=0.01):
        self.apply_diffusive_filter = apply_diffusive_filter
        self.diffusive_filter_radius = diffusive_filter_radius
        if (self.apply_diffusive_filter):
            self._BuildNodesConnectivity()

    def _ComputeDesignParameterProjectiveFilterUtilities(self):
        self._InitializeProjectiveFilter()
        self.design_parameter_projected = np.zeros(self.n_nodes)
        self.design_parameter_projected_derivatives = np.ones(self.n_nodes)

    def _InitializeProjectiveFilter(self, apply_projective_filter=True, min_mean = 0.2, max_mean = 0.5, min_projection_slope = 1e-10, max_projection_slope = 2.0, projection_change=1e-3):
        self.apply_projective_filter = apply_projective_filter
        self.projective_filter_min_mean = min_mean
        self.projective_filter_max_mean = max_mean
        self.projective_filter_min_projection_slope = min_projection_slope
        self.projective_filter_max_projection_slope = max_projection_slope
        self.projective_filter_activation_change = projection_change
        self.projective_filter_slope = 1e-15

    def _BuildNodesConnectivity(self):
        if (self.diffusive_filter_radius < 1e-10): #ensures that if no filter is imposed, at least the node itself is in neighboring nodes
            self.diffusive_filter_radius = 1e-10
        print("--|" + self.topology_optimization_stage_str + "| ---> Build Nodes Connectivity Matrix")
        mp = self._GetComputingModelPart()
        points = self._GetNodesSetCoordinates(mp.Nodes)
        nodes_tree = KDTree(points)
        self.nodes_connectivity_matrix = nodes_tree.sparse_distance_matrix(nodes_tree, self.diffusive_filter_radius, output_type="dok_matrix").tocsr()
        self.nodes_connectivity_matrix.data *= -1
        self.nodes_connectivity_matrix.data += self.diffusive_filter_radius
        self.nodes_connectivity_weigths_sum = np.array(self.nodes_connectivity_matrix.sum(axis=1)).flatten()  # Sum of each row as a 1D array
        # Normalization step: Divide non-zero entries by the corresponding row sum
        row_indices = np.repeat(np.arange(self.nodes_connectivity_matrix.shape[0]), np.diff(self.nodes_connectivity_matrix.indptr))
        self.nodes_connectivity_matrix.data /= self.nodes_connectivity_weigths_sum[row_indices]

    def _ApplyDesignParameterDiffusiveFilter(self, design_parameter):
        print("--|" + self.topology_optimization_stage_str + "| --> Apply Diffusive Filter:", self.apply_diffusive_filter)
        self.design_parameter_filtered = design_parameter
        mask = self._GetOptimizationDomainNodesMask()
        self.design_parameter_filtered[mask] = self._ApplyDiffusiveFilter(design_parameter)[mask]

    def _ApplyDesignParameterProjectiveFilter(self, design_parameter):
        print("--|" + self.topology_optimization_stage_str + "| --> Apply Projective Filter:", self.apply_projective_filter)
        mask = self._GetOptimizationDomainNodesMask()
        self.design_parameter_projected = design_parameter
        self.design_parameter_projected_derivatives = np.ones(self.n_nodes)
        if ((self.opt_it > 10) and (self.apply_projective_filter)):
            design_change_ratio = 1.0-min(1.0, max(0.0, self.design_parameter_change-self.design_parameter_change_toll)/(self.projective_filter_activation_change-self.design_parameter_change_toll))
            if (abs(design_change_ratio) < 1e-15): # design_change_ratio == 0
                new_projection_slope = self.projective_filter_min_projection_slope
            else:
                print("--|" + self.topology_optimization_stage_str + "| --> Apply Projective Filter: ---> APPLY PROJECTION")
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
                
    def _ApplyDiffusiveFilter(self, scalar_variable):
        if (self.apply_diffusive_filter):
            return (self.nodes_connectivity_matrix @ scalar_variable)
        else:
            return scalar_variable
    
    def _ApplyDiffusiveFilterDerivative(self, scalar_variable_derivative):
        if (self.apply_diffusive_filter):
            return (scalar_variable_derivative.T @ self.nodes_connectivity_matrix).T
        else:
            return scalar_variable_derivative
        

    def _UpdateOptimizationProblemPhysics(self, optimization_domain_design_parameter):
        self.design_parameter = self._InsertDesignParameterFromOptimizationDomain(optimization_domain_design_parameter)
        self._SolveTopologyOptimizationStepPhysics()
        self._EvaluateOptimizationProblem(print_results=True)
        return self.functional, self._ExtractVariableInOptimizationDomain(self.functional_derivatives_wrt_design), self.constraints, self.constraints_derivatives_wrt_design
        
    def _GetDesignParameterVariable(self):
        mp = self._GetComputingModelPart()
        design_parameter = np.zeros(self.n_nodes)
        count = 0
        for node in mp.Nodes:
            design_parameter[count] = node.GetValue(KratosMultiphysics.DESIGN_PARAMETER)
            count += 1
        return design_parameter
    
    def OutputSolutionStep(self):
        print("\n--| PRINT SOLUTION STEP OUTPUT TO FILES")
        super().OutputSolutionStep()

    def _PrintSolution(self):
        self.OutputSolutionStep()
        self._PrintFunctionalsToFile()

    def _InitializeRemeshing(self, level_set = 0.15, enable_remeshing_process = False, min_size = 0.001, max_size = 1.0):
        self.remeshing_levelset = level_set
        self.enable_remeshing = enable_remeshing_process
        if (self.IsRemeshingEnabled()):
            print("--|" + self.topology_optimization_stage_str + "| INITIALIZE REMESHING PROCESS")
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
                                                    "minimal_size"                      : """ + str(min_size) + """,
                                                    "maximal_size"                      : """ + str(max_size) + """,
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
            print("--|" + self.topology_optimization_stage_str + "| DOMAIN REMESHING")
            self.find_nodal_h.Execute()
            self.local_hessian.Execute()
            self.mmg_process.Execute()
            self._PreprocessRemeshedGeometry()

    def _PreprocessRemeshedGeometry(self):
        # GEOMETRICAL PREPROCESSING
        self._OrderNodes()
        self._OrderElements()
        self._ComputeNodalDomainSizes()

        # OPTIMIZATION PREPROCESSING
        self._OptimizationGeometricalPreprocessing()
        self._ResetConstraints()
        self.design_parameter = self._GetDesignParameterVariable()
        self._ResetPhysicsParameters()
        self._UpdatePhysicsParameters()
        self._PreprocessDerivatives() 

    def _GetPhysicsMainModelPartsList(self):
        return [self._GetPhysicsSolver().main_model_part]
    
    def EvaluateFunctionals(self, print_functional):
        self._EvaluateRequiredGradients()
        if (abs(self.functional_weights[0]) > 1e-10):
            self._EvaluateResistanceFunctional(print_functional)
        if (abs(self.functional_weights[1]) > 1e-10):
            self._EvaluateStrainRateFunctional(print_functional)
        if (abs(self.functional_weights[2]) > 1e-10):
            self._EvaluateVorticityFunctional(print_functional)

    def _CheckMaterialProperties(self):
        print("--|CHECK| Check Fluid Properties")
        self._GetSolver()._CheckMaterialProperties()

    def _UpdateRelevantPhysicsVariables(self):
        pass

    def _UpdateRelevantAdjointVariables(self):
        pass

    def _EvaluateRequiredGradients(self):
        pass
        