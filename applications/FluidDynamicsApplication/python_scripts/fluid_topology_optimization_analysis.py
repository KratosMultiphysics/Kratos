from sys import argv

import nlopt # import the library for the optimizers implementation in python
import numpy as np #import the numpy library

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
        This method executes the solution loop of the Fluid Topology Optimization Analysis
        """
        self.opt_it = 0
        print("\n--------------------------------------------------------")
        print(  "--| FLUID TOPOLOGY OPTIMIZATION PREPROCESSING")
        print("--------------------------------------------------------")
        self._InitializeTopologyOptimizationProblem()
        while (self.opt_it < self.max_it):
            self.opt_it = self.opt_it+1 # Advance Iteration
            print("\n--------------------------------------------------------")
            print(  "--| FLUID TOPOLOGY OPTIMIZATION SOLUTION LOOP. IT:", self.opt_it)
            print("--------------------------------------------------------")

            ## Current Optimization Step Solutions         
            print("\n--|INITIALIZE OPTIMIZATION STEP|")
            self._InitializeOptimizationStep((self.opt_it==1))
            print("\n--|NAVIER-STOKES SOLUTION|")
            self._RunStageSolutionLoop(1) # NAVIER-STOKES PROBLEM SOLUTION
            print("\n--|ADJOINT NAVIER-STOKES SOLUTION|")
            self._RunStageSolutionLoop(2) # ADJOINT NAVIER-STOKES PROBLEM SOLUTION
            
            ## Domain Design Update
            print("\n--|OPTIMIZER|")
            self._UpdateDomainDesign((self.opt_it==1))  # DESIGN PARAMETER UPDATE
            self._TestProblem()
            print("--------------------------------------------------------\n")

    def _RunStageSolutionLoop(self, problem_stage):
        """
        This method executes a single physics solution loop of the Topology Optimization problem.
        N.B.: must be called after the creation of the fluid model part
        problem_stage = 1: Navier Stokes solution
        problem_stage = 2: Adjoint Navier-Stokes solution
        """
        if ((problem_stage != 1) and (problem_stage != 2)): 
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
        Kratos.Logger.PrintInfo(self._GetSimulationName(), "Analysis -START- ")

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

    def _InitializeDomainDesignParameter(self):
        """
        This method will handle the design parameter initialization across the domain.
        """
        print("--|" + self.topology_optimization_stage_str + "| INITIALIZE DOMAIN DESIGN")
        mp = self._GetModelPart()
        for node in mp.Nodes:
            node.SetValue(KratosCFD.DESIGN_PARAMETER, 0.0)

    def _UpdateDomainDesign(self, first_iteration=False):
        """
        This method will handle the design parameter update.
        1) functional derivatives evaluation
        2) optimizer (MMA)?
        """
        self._SetTopologyOptimizationStage(3)
        print("--|" + self.topology_optimization_stage_str + "| DESIGN PARAMETER UPDATE NOT YET IMPLEMENTED!")
        self._EvaluateFunctionalValueAndDerivatives(first_iteration)
        self._SolveOptimizer()  

    def _EvaluateFunctionalValueAndDerivatives(self, first_iteration=False):
        """
        This method is used to evaluate the functional value and derivatives w.r.t. the design parameter
        """
        self._EvaluateFunctionalValue(first_iteration)
        self._EvaluateFunctionalDerivatives()
    
    def _EvaluateFunctionalValue(self, first_iteration=False):
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
        print("--|" + self.topology_optimization_stage_str + "| EVALUATE FUNCTIONAL VALUE")
        mp = self._GetModelPart()
        velocity = np.asarray(Kratos.VariableUtils().GetSolutionStepValuesVector(mp.Nodes, Kratos.VELOCITY, 0, self.dim)).reshape(self.n_nodes, self.dim)
        if (self.functional_weights[0] > 1e-10):
            self._EvaluateResistanceFunctional(velocity, first_iteration)
        if (self.functional_weights[1] > 1e-10):
            self._EvaluateStrainRateFunctional(velocity, first_iteration)
        if (self.functional_weights[2] > 1e-10):
            self._EvaluateVorticityFunctional(velocity, first_iteration)
        self.functional = np.dot(self.functional_weights, self.functionals)
        if (first_iteration):
            self.initial_functional_abs_value = self.functional
    
    def _EvaluateFunctionalDerivatives(self):
        """
        This method is used to evaluate the functional derivatives w.r.t the design parameter
        """
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

    def _SolveOptimizer(self):
        """
        This method is used to update the design parameter using an optimizer
        """
        dumb = 0
        return dumb
    
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
            # print("--|" + self.topology_optimization_stage_str + "| AVOIDED SOLUTION")
            # return False
        else:
            print("--|" + self.topology_optimization_stage_str + "| Time check outside NS or ADJ_NS solution")
            return False
        
    def PrintAnalysisStageProgressInformation(self):
        Kratos.Logger.PrintInfo(self._GetSimulationName(), "TOTAL STEP: ", self._GetModelPart().ProcessInfo[Kratos.STEP])
        if self.IsNavierStokesStage(): # NS
            Kratos.Logger.PrintInfo(self._GetSimulationName(), "NS STEP: ", self._GetModelPart().ProcessInfo[KratosCFD.FLUID_TOP_OPT_NS_STEP])
        elif self.IsAdjointNavierStokesStage(): # ADJ
            Kratos.Logger.PrintInfo(self._GetSimulationName(), "ADJ NS STEP: ", self._GetModelPart().ProcessInfo[KratosCFD.FLUID_TOP_OPT_ADJ_NS_STEP])
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
        This method executes the entire Fluid Topology Optimization Stage
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
    
    def _SetFunctionalWeights(self, weights = [1, 1, 0, 0, 0]):
        self.n_functionals = len(weights)
        self.initial_functionals_abs_values = np.zeros(self.n_functionals)
        self.functionals = np.zeros(self.n_functionals)
        self.functional_weights = self._NormalizeFunctionalWeights(weights)
        self._GetModelPart().ProcessInfo.SetValue(KratosCFD.FUNCTIONAL_WEIGHTS, self.functional_weights)

    def _NormalizeFunctionalWeights(self, weights):
        weights_sum = sum(weights)
        return np.asarray(weights)/weights_sum

    def _SetMaxIt(self, max_iteration = 1):
        self.max_it = max_iteration

    def _ReInitializePhysics(self):
        """
        This method Re-Initializes the physical values in order to be consistent with the solution of the current optimization step
        """
        print("--|" + self.topology_optimization_stage_str + "| RE-INITIALIZE PHYSICS")
        self._GetModelPart().ProcessInfo.SetValue(Kratos.TIME, 0.0)

    def _InitializeOptimizationStep(self, first_iteration=False):
        """
        This method Re-Initializes the physical values in order to be consistent with the solution of the current optimization step
        """
        self._SetTopologyOptimizationStage(0)
        if (not first_iteration):
                self._ReInitializePhysics()
        self._UpdateResistance()

    def _InitializeTopologyOptimizationProblem(self):
        """
        This method Initializes the topology optimization problem solution
        """
        self._InitializeDomainDesignParameter()
        self._TopologyOptimizationPreprocessing()
        self._PreprocessDerivatives()

    def _TopologyOptimizationPreprocessing(self):
        """
        This method preprocess all the useful values and quantities for the topology optimization solution
        """
        print("--|" + self.topology_optimization_stage_str + "| GENERAL PREPROCESSING")
        self.dim = self._GetModelPart().ProcessInfo.GetValue(Kratos.DOMAIN_SIZE)
        self.nodes_in_element = self.dim+1 # only for triangles(2D) and tetrahedra(3D)
        self._OrderNodes()
        self._OrderElements()
        self._ComputeNodalAreas()
        dumb = 0
        return dumb

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
        for elem in mp.Elements: 
            geom =  elem.GetGeometry()
            Ael = geom.Area()
            for node in geom:
                self.nodal_areas[node.Id-1] += Ael / contribution_factor

    def _PreprocessDerivatives(self):
        """
        This method preprocess the quantities for easy derivatives evaluation
        """
        print("--|" + self.topology_optimization_stage_str + "| PREPROCESS DERIVATIVES")
        self._PreprocessGradient()

    def _EvaluateResistanceFunctional(self, velocity, first_iteration=False):
        """
        This method computes the resistance functional: int_{\Omega}{\\alpha||u||^2}
        """
        print("--|" + self.topology_optimization_stage_str + "| ---> Resistance Functional")
        integrand = self.resistance * np.square(np.linalg.norm(velocity, axis=1)) #component-wise multiplication
        self.functionals[0] = np.dot(self.nodal_areas, integrand)
        if (first_iteration):
            self.initial_functionals_abs_values[0] = self.functionals[0] 
 
    def _EvaluateStrainRateFunctional(self, velocity, first_iteration=False):
        """
        This method computes the Strain-Rate functional: int_{\Omega}{\\2*mu*||1/2*[grad(u)+grad(u)^T]||^2}
        """
        print("--|" + self.topology_optimization_stage_str + "| ---> Strain-Rate Functional")
        mp = self._GetModelPart()
        elements_contributions = np.zeros(self.n_elements)
        for el in mp.Elements:
            nodal_velocity = velocity[self.element_nodes_ids[el.Id-1]]
            element_gradient = self._EvalElementalGradient(el, nodal_velocity)
            element_symmetric_gradient = 1.0/2.0*(element_gradient + element_gradient.transpose())
            element_symmetric_gradient_norm = np.linalg.norm(element_symmetric_gradient, ord='fro')
            element_viscosity = el.GetValue(Kratos.DYNAMIC_VISCOSITY)
            elements_contributions[el.Id-1] = 2*element_viscosity*(element_symmetric_gradient_norm**2) * el.GetGeometry().Area()
        self.functionals[1] = elements_contributions.sum()
        if (first_iteration):
            self.initial_functionals_abs_values[1] = self.functionals[1]

    def _EvaluateVorticityFunctional(self, velocity, first_iteration=False):
        """
        This method computes the Vorticity functional: int_{\Omega}{\\2*mu*||1/2*[grad(u)-grad(u)^T]||^2}
        """
        print("--|" + self.topology_optimization_stage_str + "| ---> Vorticity Functional")
        mp = self._GetModelPart()
        elements_contributions = np.zeros(self.n_elements)
        for el in mp.Elements:
            nodal_velocity = velocity[self.element_nodes_ids[el.Id-1]]
            element_gradient = self._EvalElementalGradient(el, nodal_velocity)
            element_antisymmetric_gradient = 1.0/2.0*(element_gradient - element_gradient.transpose())
            element_antisymmetric_gradient_norm = np.linalg.norm(element_antisymmetric_gradient, ord='fro')
            element_viscosity = el.GetValue(Kratos.DYNAMIC_VISCOSITY)
            elements_contributions[el.Id-1] = 2*element_viscosity*(element_antisymmetric_gradient_norm**2) * el.GetGeometry().Area()
        self.functionals[2] = elements_contributions.sum()
        if (first_iteration):
            self.initial_functionals_abs_values[2] = self.functionals[2]

    def _TestProblem(self):
        dumb = 0
        return dumb

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
        
    def _EvalElementalGradient(self, element, nodal_variable):
        shape_functions_derivatives = self.shape_functions_derivatives[element.Id-1]
        return (nodal_variable.transpose() @ shape_functions_derivatives) # matrix multilpication
    
    def _GetShapeFunctionsDerivatives(self, element):
        gradients = np.asarray(element.GetGeometry().ShapeFunctionDerivatives(1,0))
        jacobian  = np.asarray(element.GetGeometry().Jacobian(0))
        inv_jacobian = np.linalg.inv(jacobian)
        global_gradient = gradients @ inv_jacobian
        return global_gradient
    
    def GetObjectiveFunctionAndGradient(self, design_parameter, grad):
        functional = self.functional
        if grad.size > 0:
            grad = self.functional_derivatives_wrt_design
        return functional
    
    def EvalVolumeConstraintAndDerivative(self, design_parameter, grad):
        if grad.size > 0:
            grad[:] = 1.0 / self.nodal_areas
        return np.dot(design_parameter, self.nodal_areas)








        