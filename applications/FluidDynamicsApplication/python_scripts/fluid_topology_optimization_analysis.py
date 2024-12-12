from sys import argv

import numpy as np #import the numpy library
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
        self._SetMaxIt()  

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
        
    def _SetMaxIt(self, max_iteration = 1):
        """
        This method sets the maximum iteration for the topology optimization process
        """
        self.max_it = max_iteration

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
    
    def _InitializeResistance(self, min_value = 0.0, max_value = 1000.0, interpolation_slope = 1.0):
        print("--|" + self.topology_optimization_stage_str + "| INITIALIZE RESISTANCE")
        self.alpha_min = min_value
        self.alpha_max = max_value
        self.q = interpolation_slope
        self.resistance = np.zeros(self.n_nodes)
        self.resistance_derivative_wrt_design = np.zeros(self.n_nodes)
    
    def _UpdateResistance(self, min_value = 0.0, max_value = 1000.0, interpolation_slope = 1.0):
        """
        This method will handle the resistance update, that will be based on the value of the design parameter.
        Now it is used to "play" with the resistance
        """
        print("--|" + self.topology_optimization_stage_str + "| UPDATE RESISTANCE")
        if (not hasattr(self, 'alpha_min')) or (not hasattr(self, 'alpha_max') or (not hasattr(self, 'q'))):
            self._InitializeResistance(min_value, max_value, interpolation_slope)
        mp = self._GetModelPart()
        for node in mp.Nodes:
                nodal_gamma = node.GetValue(KratosCFD.DESIGN_PARAMETER)
                nodal_resistance = self.alpha_min + (self.alpha_max-self.alpha_min)*(self.q*nodal_gamma)/(self.q+1-nodal_gamma)
                node.SetValue(KratosCFD.RESISTANCE, nodal_resistance)
                self.resistance[node.Id-1] = nodal_resistance
                self.resistance_derivative_wrt_design[node.Id-1] = (self.alpha_max-self.alpha_min)*(self.q*(self.q+1))/((self.q+1-nodal_gamma)*(self.q+1-nodal_gamma))
       
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
        if self._GetModelPart().ProcessInfo[Kratos.IS_RESTARTED]:
            self.time = self._GetModelPart().ProcessInfo[Kratos.TIME]
        else:
            self.time = self.project_parameters["problem_data"]["start_time"].GetDouble()
            self._GetModelPart().ProcessInfo[Kratos.TIME] = self.time
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
        self.initial_functionals_abs_values = np.zeros(self.n_functionals)
        self.functionals = np.zeros(self.n_functionals)
        self.functional_weights = self._NormalizeFunctionalWeights(np.asarray(weights))
        self._GetModelPart().ProcessInfo.SetValue(KratosCFD.FUNCTIONAL_WEIGHTS, self.functional_weights)

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
        self._GeometricalPreprocessing()
        self._InitializeOptimization()
        self._InitializeDomainDesign()
        self._PreprocessDerivatives() 

    def SolveFluidTopologyOptimization(self):
        while (self.opt_it < self.max_it):
            self.opt_it = self.opt_it+1
            print("\n--------------------------------------------------------")
            print(  "--| FLUID TOPOLOGY OPTIMIZATION SOLUTION LOOP. IT:", self.opt_it)
            print("--------------------------------------------------------")
            self._SolveOptimizer(self.opt_it == 1)
            self._SolveTopologyOptimizationStepPhysics()
            self._EvaluateOptimizationProblem((self.opt_it==1))

    def _SolveOptimizer(self, first_iteration = False):
        print("\n--|OPTIMIZER|")
        self._SetTopologyOptimizationStage(3)
        if (first_iteration):
            print("--|" + self.topology_optimization_stage_str + "| First Iteration: DO NOTHING")
            pass
        else:
            print("--|" + self.topology_optimization_stage_str + "| UPDATE DESIGN PARAMETER")
            self._SolveMMA(self.design_parameter, self.n_design_parameters, self.n_optimization_constraints)
    
    def _SolveTopologyOptimizationStepPhysics(self):      
        ## Current Optimization Step Solutions 
        self._InitializeTopologyOptimizationStepPhysicsSolution((self.opt_it==1))
        self._SolveNaviersStokesProblem() # NAVIER-STOKES PROBLEM SOLUTION
        self._SolveAdjointNaviersStokesProblem() # ADJOINT NAVIER-STOKES PROBLEM SOLUTION   
    
    def _GeometricalPreprocessing(self):
        """
        This method preprocess all the useful values and quantities for the topology optimization solution
        """
        print("--|" + self.topology_optimization_stage_str + "| GEOMETRICAL PREPROCESSING")
        mp = self._GetModelPart()
        self.dim = mp.ProcessInfo.GetValue(Kratos.DOMAIN_SIZE)
        self.nodes_in_element = self.dim+1 # only for triangles(2D) and tetrahedra(3D)
        self._OrderNodes()
        self._OrderElements()
        self._ComputeNodalAreas()
    
    def _InitializeOptimization(self):  
        print("--|" + self.topology_optimization_stage_str + "| INITIALIZE OPTIMIZATION")
        self.opt_it = 0
        self.n_design_parameters = self.n_nodes
        self.n_optimization_constraints = 1  
        self.constraints = np.zeros(self.n_optimization_constraints)  
        self.constraints_derivatives_wrt_design = np.zeros((self.n_optimization_constraints, self.n_design_parameters))

    def _InitializeDomainDesign(self):
        print("--|" + self.topology_optimization_stage_str + "| INITIALIZE DOMAIN DESIGN")
        self._SetMaxDomainVolumeFraction()
        self._InitializeDomainDesignParameter()

    def _SetMaxDomainVolumeFraction(self, volume_fraction = 1.0):
        self.max_volume_fraction = volume_fraction

    def _InitializeDomainDesignParameter(self):
        """
        This method will handle the design parameter initialization across the domain.
        """
        self.design_parameter = np.zeros(self.n_nodes)
        mp = self._GetModelPart()
        for node in mp.Nodes:
            node.SetValue(KratosCFD.DESIGN_PARAMETER, self.design_parameter[node.Id-1])

    def _OrderNodes(self):
        """
        This method orders the model part nodes in increasing order starting from 1
        """
        count = 0
        mp = self._GetModelPart()
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
        count = 0
        mp = self._GetModelPart()
        for el in mp.Elements:
            count = count+1
            el.Id = count
        if (count != len(mp.Elements)):
            print("\nERROR: wrong reordering of elements ids. The counted number of elements is different from len(mp.Elements)\n")
        self.n_elements = count

    def _ComputeNodalAreas(self):
        """
        This method compute the nodal area - not vectorized but it is done once. WORKS ONLY FOR TRIANGULAR AND TETRAHEDRAL MESH
        """
        mp = self._GetModelPart()
        contribution_factor = 1.0/self.nodes_in_element 
        self.nodal_areas = np.zeros(len(mp.Nodes))
        self.elemental_area = np.zeros(len(mp.Elements))
        self.total_area = 0.0
        for elem in mp.Elements: 
            geom =  elem.GetGeometry()
            area_elem = geom.Area()
            self.elemental_area[elem.Id-1] = area_elem
            self.total_area = self.total_area + area_elem
            for node in geom:
                self.nodal_areas[node.Id-1] += area_elem * contribution_factor

    def _PreprocessDerivatives(self):
        """
        This method preprocess the quantities for easy derivatives evaluation
        """
        print("--|" + self.topology_optimization_stage_str + "| PREPROCESS DERIVATIVES")
        self._PreprocessGradient()

    def _PreprocessGradient(self):
        self.shape_functions_derivatives = np.zeros((self.n_elements, self.nodes_in_element, self.dim))
        self.element_nodes_ids = np.zeros((self.n_elements, self.nodes_in_element), dtype=int)
        mp = self._GetModelPart()
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
        self.design_parameter = design_parameter
        self._UpdateDesignParameterVariable()

    def _UpdateDesignParameterVariable(self):
        mp = self._GetModelPart()
        for node in mp.Nodes:
                node.SetValue(KratosCFD.DESIGN_PARAMETER, self.design_parameter[node.Id-1])
    
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
            self.OutputSolutionStep()
        print("--|" + top_opt_stage_str + "| END SOLUTION LOOP")

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
        self._GetModelPart().ProcessInfo.SetValue(Kratos.TIME, 0.0)

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
        move = 1.0
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
        self._UpdateDesignParameterAndResistance(xval.flatten())

    def _EvaluateOptimizationProblem(self, first_iteration = False):
        print("\n--|EVALUATE OPTIMIZATION PROBLEM|")
        self._UpdateOptimizationProblem(self.design_parameter, first_iteration)
        self._PrintOptimizationProblem()
    
    def _UpdateOptimizationProblem(self, design_parameter, first_iteration = False):
        self._EvaluateFunctionalAndDerivatives(design_parameter, first_iteration)
        self._EvaluateConstraintsAndDerivatives(design_parameter)
        return self.functional, self.functional_derivatives_wrt_design, self.constraints, self.constraints_derivatives_wrt_design
    
    def _EvaluateConstraintsAndDerivatives(self, design_parameter):
         self._EvaluateVolumeConstraintAndDerivative(design_parameter)
         self.constraints[0] = self.volume_constraint
         self.constraints_derivatives_wrt_design[0,:] = self.volume_constraint_derivatives_wrt_design

    def _EvaluateFunctionalAndDerivatives(self, design_parameter, first_iteration=False):
        """
        This method is used to evaluate the functional value and derivatives w.r.t. the design parameter
        """
        self._SetTopologyOptimizationStage(3)
        self._UpdateDesignParameterAndResistance(design_parameter)
        self._EvaluateFunctional(first_iteration)
        self._EvaluateFunctionalDerivatives()
    
    def _EvaluateFunctional(self, first_iteration=False):
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
        mp = self._GetModelPart()
        velocity = np.asarray(Kratos.VariableUtils().GetSolutionStepValuesVector(mp.Nodes, Kratos.VELOCITY, 0, self.dim)).reshape(self.n_nodes, self.dim)
        if (abs(self.functional_weights[0]) > 1e-10):
            self._EvaluateResistanceFunctional(velocity, first_iteration)
        if (abs(self.functional_weights[1]) > 1e-10):
            self._EvaluateStrainRateFunctional(velocity, first_iteration)
        if (abs(self.functional_weights[2]) > 1e-10):
            self._EvaluateVorticityFunctional(velocity, first_iteration)
        self.functional = np.dot(self.functional_weights, self.functionals)
        self.weighted_functionals = self.functional_weights * self.functionals
        if (first_iteration):
            self.initial_functional_abs_value = self.functional
            self.initial_functionals_abs_value = self.functionals
            self.initial_weighted_functionals_abs_value = self.weighted_functionals

    def _EvaluateResistanceFunctional(self, velocity, first_iteration=False, print_functional=False):
        """
        This method computes the resistance functional: int_{\Omega}{\\alpha||u||^2}
        """
        print("--|" + self.topology_optimization_stage_str + "| --> Resistance Functional")
        nodal_velocity_norm = np.linalg.norm(velocity, axis=1)
        integrand = self.resistance * np.square(nodal_velocity_norm) #component-wise multiplication
        self.functionals[0] = np.dot(self.nodal_areas, integrand)
        if (first_iteration):
            self.initial_functionals_abs_values[0] = self.functionals[0] 
        if (print_functional):
            print("--|" + self.topology_optimization_stage_str + "| ---> Resistance Functional (no weight):", self.functionals[0])
 
    def _EvaluateStrainRateFunctional(self, velocity, first_iteration=False, print_functional=False):
        """
        This method computes the Strain-Rate functional: int_{\Omega}{\\2*mu*||1/2*[grad(u)+grad(u)^T]||^2}
        """
        print("--|" + self.topology_optimization_stage_str + "| --> Strain-Rate Functional")
        v = velocity[self.element_nodes_ids[:]]
        g = np.matmul(np.transpose(v, axes=(0,2,1)), self.shape_functions_derivatives)
        g_sym = 1.0/2.0 * (g+(np.transpose(g, axes=(0,2,1))))
        g_sym_norm_sq = (np.linalg.norm(g_sym, ord='fro', axis=(1, 2)))**2.0
        mu = self._GetModelPart().Elements[1].Properties.GetValue(Kratos.DYNAMIC_VISCOSITY)
        self.functionals[1] = 2*mu* np.dot(g_sym_norm_sq, self.elemental_area)
        if (first_iteration):
            self.initial_functionals_abs_values[1] = self.functionals[1]
        if (print_functional):
            print("--|" + self.topology_optimization_stage_str + "| ---> Strain-Rate Functional (no weight):", self.functionals[1])
    
    def _EvaluateVorticityFunctional(self, velocity, first_iteration=False, print_functional=False):
        """
        This method computes the Vorticity functional: int_{\Omega}{\\2*mu*||1/2*[grad(u)-grad(u)^T]||^2}
        """
        print("--|" + self.topology_optimization_stage_str + "| --> Vorticity Functional")
        v = velocity[self.element_nodes_ids[:]]
        g = np.matmul(np.transpose(v, axes=(0,2,1)), self.shape_functions_derivatives)
        g_sym = 1.0/2.0 * (g+(np.transpose(g, axes=(0,2,1))))
        g_sym_norm_sq = (np.linalg.norm(g_sym, ord='fro', axis=(1, 2)))**2.0
        mu = self._GetModelPart().Elements[1].Properties.GetValue(Kratos.DYNAMIC_VISCOSITY)
        self.functionals[1] = 2*mu* np.dot(g_sym_norm_sq, self.elemental_area)
        if (first_iteration):
            self.initial_functionals_abs_values[2] = self.functionals[2]
        if (print_functional):
            print("--|" + self.topology_optimization_stage_str + "| ---> Vorticity Functional: (no weight)", self.functionals[2])    
    
    def _EvaluateFunctionalDerivatives(self, first_it = False):
        """
        This method is used to evaluate the functional derivatives w.r.t the design parameter
        """
        self._SetTopologyOptimizationStage(3)
        print("--|" + self.topology_optimization_stage_str + "| EVALUATE FUNCTIONAL DERIVATIVES")
        mp = self._GetModelPart()
        velocity = np.asarray(Kratos.VariableUtils().GetSolutionStepValuesVector(mp.Nodes, Kratos.VELOCITY, 0, self.dim)).reshape(self.n_nodes, self.dim)
        velocity_adjoint= np.asarray(Kratos.VariableUtils().GetSolutionStepValuesVector(mp.Nodes, KratosCFD.VELOCITY_ADJ, 0, self.dim)).reshape(self.n_nodes, self.dim)
        # NEXT LINE EXPLANATION:
        # df/dgamma_i = f_weight[0]*dAlpha/dgamma[i]*||vel[i]||^2 + dAlpha/dgamma[i]*(vel[i].dot(vel_adj[i])) = dAlpha/dgamma[i]*vel[i].dot(f_weight[0]*vel[i] + vel_adj[i])
        # then component by component is multiplied by the nodal area of influence
        # !!!!!!!! ASK RICCARDO IF IT IS CORRECT !!!!!!!!!
        self.functional_derivatives_wrt_design = self.resistance_derivative_wrt_design * np.sum(velocity * (self.functional_weights[0]*velocity + velocity_adjoint), axis=1) * self.nodal_areas
        for node in mp.Nodes:
            node.SetValue(KratosCFD.FUNCTIONAL_DERIVATIVE, self.functional_derivatives_wrt_design[node.Id-1])
        if (first_it):
            self.initial_functional_derivatives_wrt_design = self.functional_derivatives_wrt_design
        self.functional_derivatives_wrt_design = np.asarray([self.functional_derivatives_wrt_design]).T

    def _EvaluateVolumeConstraintAndDerivative(self, design_parameter):
        self.volume_fraction = 1.0 - np.dot(design_parameter, self.nodal_areas)/self.total_area
        self.volume_constraint = self.volume_fraction - self.max_volume_fraction
        self.volume_constraint_derivatives_wrt_design = -1.0 * self.nodal_areas / self.total_area
    
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
        self._GetModelPart().ProcessInfo.SetValue(KratosCFD.FLUID_TOP_OPT_PROBLEM_STAGE, self.topology_optimization_stage)

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
    
    def _GetModelPart(self):
        """
        This method returns the computing model part for the current solver
        """    
        return self._GetSolver().GetComputingModelPart()

    def _GetNodes(self):
        """
        This method returns the nodes for the current solver
        """
        return self._GetModelPart().Nodes
    
    def _GetElements(self):
        """
        This method returns the elements for the current solver
        """
        return self._GetModelPart().Elements
    
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
    
    ## PRINTS
    def _PrintOptimizationProblem(self):
        print("\n--|PRINT OPTIMIZATION PROBLEM DATA|")
        self._PrintFunctionals()
        self._PrintConstraints()
    
    def _PrintFunctionals(self):
        print("--|" + self.topology_optimization_stage_str + "| TOTAL FUNCTIONAL:", self.functional)
        if (abs(self.functional_weights[0]) > 1e-10):
            print("--|" + self.topology_optimization_stage_str + "| ---> Resistance Functional (" + str(self.functional_weights[0]) + "):", self.weighted_functionals[0])
        if (abs(self.functional_weights[1]) > 1e-10):
            print("--|" + self.topology_optimization_stage_str + "| ---> Strain-Rate Functional (" + str(self.functional_weights[1]) + "):", self.weighted_functionals[1])
        if (abs(self.functional_weights[2]) > 1e-10):
            print("--|" + self.topology_optimization_stage_str + "| ---> Vorticity Functional (" + str(self.functional_weights[2]) + "):", self.weighted_functionals[2])
        
    def _PrintConstraints(self):
        self._PrintVolumeConstraint()
    
    def _PrintVolumeConstraint(self):
        print("--|" + self.topology_optimization_stage_str + "| VOLUME CONSTRAINT:", self.volume_constraint)
        print("--|" + self.topology_optimization_stage_str + "| ---> Volume Fraction:", self.volume_fraction)

    def PrintAnalysisStageProgressInformation(self):
        Kratos.Logger.PrintInfo(self._GetSimulationName(), "TOTAL STEP: ", self._GetModelPart().ProcessInfo[Kratos.STEP])
        if self.IsNavierStokesStage(): # NS
            Kratos.Logger.PrintInfo(self._GetSimulationName(), "NS STEP: ", self._GetModelPart().ProcessInfo[KratosCFD.FLUID_TOP_OPT_NS_STEP])
        elif self.IsAdjointNavierStokesStage(): # ADJ
            Kratos.Logger.PrintInfo(self._GetSimulationName(), "ADJ NS STEP: ", self._GetModelPart().ProcessInfo[KratosCFD.FLUID_TOP_OPT_ADJ_NS_STEP])
        else:
            Kratos.Logger.PrintInfo(self._GetSimulationName(), "Invalid Value of the Topology Optimization Stage. TOP_OPT_STAGE: ", self.topology_optimization_stage)
        Kratos.Logger.PrintInfo(self._GetSimulationName(), "TIME: ", self.time)













        