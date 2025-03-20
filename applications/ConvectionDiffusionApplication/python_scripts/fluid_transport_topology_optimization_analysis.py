from sys import argv

import numpy as np #import the numpy library
import scipy as sp #import the scipy library
from scipy.spatial import KDTree
from scipy.sparse import dok_matrix, lil_matrix
import mmapy as MMA # import the MMA subroutines python library: https://github.com/arjendeetman/GCMMA-MMA-Python
import time

import KratosMultiphysics as KratosMultiphysics
# Import applications
import KratosMultiphysics.FluidDynamicsApplication as KratosCFD
import KratosMultiphysics.ConvectionDiffusionApplication as KratosCD
import KratosMultiphysics.MeshingApplication as KratosMMG

from KratosMultiphysics.FluidDynamicsApplication.fluid_dynamics_analysis import FluidDynamicsAnalysis
from KratosMultiphysics.ConvectionDiffusionApplication.convection_diffusion_analysis import ConvectionDiffusionAnalysis
from KratosMultiphysics.FluidDynamicsApplication.fluid_topology_optimization_analysis import FluidTopologyOptimizationAnalysis
from KratosMultiphysics.ConvectionDiffusionApplication.transport_topology_optimization_analysis import TransportTopologyOptimizationAnalysis
from KratosMultiphysics.analysis_stage import AnalysisStage
from KratosMultiphysics.FluidDynamicsApplication import fluid_topology_optimization_solver
from KratosMultiphysics.ConvectionDiffusionApplication import transport_topology_optimization_solver
from KratosMultiphysics.ConvectionDiffusionApplication import fluid_transport_topology_optimization_solver

class FluidTransportTopologyOptimizationAnalysis(TransportTopologyOptimizationAnalysis):
    def __init__(self,model,parameters):
        super().__init__(model,parameters) 

    def _SetTopologyOptimizationName(self):
        self.topology_optimization_name = "FLUID-TRANSPORT"

    def _CreateTopologyOptimizationSolvers(self):
        """
        This method creates the NS and ADJ_NS solvers
        """
        self.physics_solver = fluid_transport_topology_optimization_solver.CreateSolver(self.model, self.project_parameters)
        self.adjoint_solver = fluid_transport_topology_optimization_solver.CreateSolver(self.model, self.project_parameters, isAdjointSolver=True)

    def _CreateSolver(self, isAdjointSolver = False):
        """
        This method creates a solver
        isAdjointSolver == False --> physics_solver
        isAdjointSolver == True  --> adjoint_solver
        """
        return fluid_transport_topology_optimization_solver.CreateSolver(self.model, self.project_parameters, isAdjointSolver)
    
    def PrepareSolvers(self):
        """This method prepares the NS and ADJ_NS Solvers in the AnalysisStage 
        Usage: It is designed to be called ONCE, BEFORE the execution of the solution-loop
        Prepare Solvers : ImportModelPart -> PrepareModelPart -> AddDofs
        """
        self.PreparePhysicsSolver()
        self.PrepareAdjointSolver()
    
    def _InitializeTopologyOptimizationProblem(self):
        super()._InitializeTopologyOptimizationProblem()

    def _RunStageSolutionLoop(self):
        """
        This method executes a single physics solution of the Topology Optimization problem loop
        N.B.: must be called after the creation of the fluid-transport model part
        problem_stage = 1: Navier Stokes solution
        problem_stage = 2: Adjoint Navier-Stokes solution
        """
        self._PrepareTransportSettings()
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

    def _PrepareTransportSettings(self):
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
            print("\n ERROR: _PrepareTransportSettings in the wrong topology optimization stage.\n")
        
    def _GetPhysicsMainModelPartsList(self):
        return [self._GetFluidSolver().main_model_part, self._GetTransportSolver().main_model_part]

    def PrintAnalysisStageProgressInformation(self):
        KratosMultiphysics.Logger.PrintInfo(self._GetSimulationName(), "TOTAL STEP: ", self._GetComputingModelPart().ProcessInfo[KratosMultiphysics.STEP])
        if self.IsPhysicsStage(): # NS
            KratosMultiphysics.Logger.PrintInfo(self._GetSimulationName(), "NS STEP: ", self._GetComputingModelPart().ProcessInfo[KratosCFD.FLUID_TOP_OPT_NS_STEP])
            KratosMultiphysics.Logger.PrintInfo(self._GetSimulationName(), " T STEP: ", self._GetComputingModelPart().ProcessInfo[KratosCD.TRANSPORT_TOP_OPT_T_STEP])
        elif self.IsAdjointStage(): # ADJ
            KratosMultiphysics.Logger.PrintInfo(self._GetSimulationName(), "ADJ NS STEP: ", self._GetComputingModelPart().ProcessInfo[KratosCFD.FLUID_TOP_OPT_ADJ_NS_STEP])
            KratosMultiphysics.Logger.PrintInfo(self._GetSimulationName(), "ADJ  T STEP: ", self._GetComputingModelPart().ProcessInfo[KratosCD.TRANSPORT_TOP_OPT_ADJ_T_STEP])
        else:
            KratosMultiphysics.Logger.PrintInfo(self._GetSimulationName(), "Invalid Value of the Topology Optimization Stage. TOP_OPT_STAGE: ", self.topology_optimization_stage)
        KratosMultiphysics.Logger.PrintInfo(self._GetSimulationName(), "TIME: ", self.time) 

    def _SetTransportSurfaceSourceFromFunctional(self):
        process = self._FindAdjointSurfaceSourceProcess()
        process.value += self.functional_weights[3]*self.avg_outlet_transport_scalar_diff

    def _FindAdjointSurfaceSourceProcess(self):
        # Get the list of processes
        process_list = self._GetListOfProcesses()
        # Find a process by class name
        for process in process_list:
            if isinstance(process, KratosMultiphysics.assign_scalar_variable_process.AssignScalarVariableProcess):
                if (process.variable == KratosCD.FACE_HEAT_FLUX_ADJ):
                    found = True
                    surface_source_process = process
        if (found):
            return surface_source_process
        else:
            KratosMultiphysics.Logger.PrintError("Not Found Adjoint Volume Source Process", "It should be using the 'FACE_HEAT_FLUX_ADJ' variable")
    
    def _FindAdjointVolumeSourceProcess(self):
        # Get the list of processes
        process_list = self._GetListOfProcesses()
        # Find a process by class name
        found = False
        for process in process_list:
            if isinstance(process, KratosMultiphysics.assign_scalar_variable_process.AssignScalarVariableProcess):
                if (process.variable == KratosCD.HEAT_FLUX_ADJ):
                    found = True
                    volume_source_process = process
        if (found):
            return volume_source_process
        else:
            KratosMultiphysics.Logger.PrintError("Not Found Adjoint Volume Source Process", "It should be using the 'HEAT_FLUX_ADJ' variable")

    def _PrintFunctionals(self):
        print("--|" + self.topology_optimization_stage_str + "| TOTAL FUNCTIONAL  :", self.functional)
        print("--|" + self.topology_optimization_stage_str + "| INITIAL FUNCTIONAL:", self.initial_functional)
        if (abs(self.functional_weights[0]) > 1e-10):
            print("--|" + self.topology_optimization_stage_str + "| ---> Resistance Functional (" + str(self.functional_weights[0]) + "):", self.weighted_functionals[0])
        if (abs(self.functional_weights[1]) > 1e-10):
            print("--|" + self.topology_optimization_stage_str + "| ---> Strain-Rate Functional (" + str(self.functional_weights[1]) + "):", self.weighted_functionals[1])
        if (abs(self.functional_weights[2]) > 1e-10):
            print("--|" + self.topology_optimization_stage_str + "| ---> Vorticity Functional (" + str(self.functional_weights[2]) + "):", self.weighted_functionals[2])
        if (abs(self.functional_weights[3]) > 1e-10):
            print("--|" + self.topology_optimization_stage_str + "| ---> Outlet Transport Scalar Functional (" + str(self.functional_weights[3]) + "):", self.weighted_functionals[3])
        if (abs(self.functional_weights[4]) > 1e-10):
            print("--|" + self.topology_optimization_stage_str + "| ---> Focus Region Transport Scalar Functional (" + str(self.functional_weights[4]) + "):", self.weighted_functionals[4])
        if (abs(self.functional_weights[5]) > 1e-10):
            print("--|" + self.topology_optimization_stage_str + "| ---> Transport Scalar Diffusion Functional (" + str(self.functional_weights[5]) + "):", self.weighted_functionals[5])
        if (abs(self.functional_weights[6]) > 1e-10):
            print("--|" + self.topology_optimization_stage_str + "| ---> Transport Scalar Convection Functional (" + str(self.functional_weights[6]) + "):", self.weighted_functionals[6])
        if (abs(self.functional_weights[7]) > 1e-10):
            print("--|" + self.topology_optimization_stage_str + "| ---> Transport Scalar Decay Functional (" + str(self.functional_weights[7]) + "):", self.weighted_functionals[7])
        if (abs(self.functional_weights[8]) > 1e-10):
            print("--|" + self.topology_optimization_stage_str + "| ---> Transport Scalar Source Functional (" + str(self.functional_weights[8]) + "):", self.weighted_functionals[8])

    def _SetFunctionalWeights(self):
        super()._SetFunctionalWeights()
        self.EvaluateTotalFunctional()
        
    def _InitializeFunctionalWeights(self):
        fluid_weights     = self._ImportFluidFunctionalWeights()
        transport_weights = self._ImportTransportFunctionalWeights()
        coupling_weights  = self._ImportCouplingFunctionalWeights()
        # normalize weights
        self.normalized_fluid_functional_weights = self._NormalizeFunctionalWeights(np.asarray(fluid_weights))
        self.normalized_transport_functional_weights = self._NormalizeFunctionalWeights(np.asarray(transport_weights))
        self.normalized_coupling_functional_weights = self._NormalizeFunctionalWeights(np.asarray(coupling_weights))
        # get number of functionals
        self.n_fluid_functionals = len(self.normalized_fluid_functional_weights)
        self.n_transport_functionals = len(self.normalized_transport_functional_weights)
        self.n_functionals = self.n_fluid_functionals+self.n_transport_functionals
        # initialize initial functionals vector container
        self.initial_fluid_functionals_values = np.zeros(self.n_fluid_functionals)
        self.initial_transport_functionals_values = np.zeros(self.n_transport_functionals)
        self.initial_functionals_values = np.zeros(self.n_functionals)
        # initialize functionals vector container
        self.fluid_functionals = np.zeros(self.n_fluid_functionals)
        self.transport_functionals = np.zeros(self.n_transport_functionals)
        self.functionals = np.zeros(self.n_functionals)
        self.EvaluateFunctionals(print_functional=False)
        self._SetInitialFunctionals()
        self.functional_weights = self._RescaleFunctionalWeightsByInitialValues()

    def _ImportCouplingFunctionalWeights(self):
        coupling_weights = [0.0, 0.0]
        functional_weights_parameters = self.optimization_parameters["optimization_settings"]["optimization_problem_settings"]["functional_weights"]["coupling"]
        coupling_weights[0] = functional_weights_parameters["fluid"].GetDouble()
        coupling_weights[1] = functional_weights_parameters["transport"].GetDouble()
        return coupling_weights

    def _SetInitialFunctionals(self):
        self.initial_fluid_functional = np.dot(self.normalized_fluid_functional_weights, self.initial_fluid_functionals_values)
        self.initial_transport_functional = np.dot(self.normalized_transport_functional_weights, self.initial_transport_functionals_values)
        if (abs(self.initial_fluid_functional) < 1e-10):
            print("[WARNING] Initial fluid functional is zero")
        if (abs(self.initial_transport_functional) < 1e-10):
            print("[WARNING] Initial transport functional is zero")
        self.initial_coupling_functionals = np.asarray([self.initial_fluid_functional, self.initial_transport_functional])
        self.initial_coupling_functionals_abs_value = np.abs(self.initial_coupling_functionals)

    def _RescaleFunctionalWeightsByInitialValues(self):
        # fluid
        if ((np.abs(self.normalized_coupling_functional_weights[0]) < 1e-10) or (np.sum(np.abs(self.normalized_fluid_functional_weights)) < 1e-10)):
            fluid_functional_weights = np.zeros(self.normalized_fluid_functional_weights.size)
        else:
            fluid_functional_weights  = self.normalized_fluid_functional_weights * self.normalized_coupling_functional_weights[0]
            if (abs(self.initial_fluid_functional) > 1e-15):
                fluid_functional_weights /= abs(self.initial_fluid_functional)
        # transport
        if ((np.abs(self.normalized_coupling_functional_weights[1]) < 1e-10) or (np.sum(np.abs(self.normalized_transport_functional_weights)) < 1e-10)):
            transport_functional_weights = np.zeros(self.normalized_transport_functional_weights.size)
        else:
            transport_functional_weights  = self.normalized_transport_functional_weights * self.normalized_coupling_functional_weights[1]
            if (abs(self.initial_transport_functional) > 1e-15):
                transport_functional_weights /= abs(self.initial_transport_functional)
        return np.concatenate((fluid_functional_weights, transport_functional_weights))

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
        self.EvaluateTotalFunctional()

    def _EvaluateRequiredGradients(self):
        super()._EvaluateRequiredGradients()
        self._ComputeScalarVariableNodalGradient(KratosMultiphysics.TEMPERATURE, KratosMultiphysics.TEMPERATURE_GRADIENT)
        self._ComputeScalarVariableNodalGradient(KratosMultiphysics.TEMPERATURE_ADJ, KratosMultiphysics.TEMPERATURE_ADJ_GRADIENT)

    def EvaluateFunctionals(self, print_functional):
        self._EvaluateRequiredGradients()
        # FLUID FUNCTIONALS
        if (abs(self.normalized_fluid_functional_weights[0]) > 1e-10):
            self._EvaluateResistanceFunctional(print_functional)
        if (abs(self.normalized_fluid_functional_weights[1]) > 1e-10):
            self._EvaluateStrainRateFunctional(print_functional)
        if (abs(self.normalized_fluid_functional_weights[2]) > 1e-10):
            self._EvaluateVorticityFunctional(print_functional)
        # TRANSPORT FUNCTIONALS
        if (abs(self.normalized_transport_functional_weights[0]) > 1e-10):
            self._EvaluateOutletTransportScalarFunctional(print_functional)
        if (abs(self.normalized_transport_functional_weights[1]) > 1e-10):
            self._EvaluateFocusRegionTransportScalarFunctional(print_functional)
        if (abs(self.normalized_transport_functional_weights[2]) > 1e-10):
            self._EvaluateTransportScalarDiffusionFunctional(print_functional)
        if (abs(self.normalized_transport_functional_weights[3]) > 1e-10):
            self._EvaluateTransportScalarConvectionFunctional(print_functional)
        if (abs(self.normalized_transport_functional_weights[4]) > 1e-10):
            self._EvaluateTransportScalarDecayFunctional(print_functional)
        if (abs(self.normalized_transport_functional_weights[5]) > 1e-10):
            self._EvaluateTransportScalarSourceFunctional(print_functional)
            

    def EvaluateTotalFunctional(self):
        self.functionals = np.concatenate((self.fluid_functionals, self.transport_functionals))
        self.weighted_functionals = self.functional_weights * self.functionals
        self.functional = np.sum(self.weighted_functionals)
        if (self.first_iteration):
            self.initial_functional = self.functional

    def _EvaluateResistanceFunctional(self, print_functional=False):
        """
        This method computes the resistance functional: int_{\Omega}{\\alpha||u||^2}
        """
        mp = self._GetComputingModelPart()
        velocity = np.asarray(KratosMultiphysics.VariableUtils().GetSolutionStepValuesVector(mp.Nodes, KratosMultiphysics.VELOCITY, 0, self.dim)).reshape(self.n_nodes, self.dim)
        nodal_velocity_norm = np.linalg.norm(velocity, axis=1)
        integrand = self.resistance * (nodal_velocity_norm**2) #component-wise multiplication
        self.fluid_functionals[0] = np.dot(self.nodal_domain_sizes, integrand)
        if (self.first_iteration):
            self.initial_fluid_functionals_values[0] = self.fluid_functionals[0] 
        if (print_functional):
            print("--|" + self.topology_optimization_stage_str + "| ---> Resistance Functional (no weight):", self.fluid_functionals[0])
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
        self.fluid_functionals[1] = 2*mu* np.dot(vel_symmetric_gradient_norm_squared, self.elemental_domain_size)
        if (self.first_iteration):
            self.initial_fluid_functionals_values[1] = self.fluid_functionals[1]
        if (print_functional):
            print("--|" + self.topology_optimization_stage_str + "| ---> Strain-Rate Functional (no weight):", self.fluid_functionals[1])
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
        self.fluid_functionals[2] = 2*mu* np.dot(vel_antisymmetric_gradient_norm_squared, self.elemental_domain_size)
        if (self.first_iteration):
            self.initial_fluid_functionals_values[2] = self.fluid_functionals[2]
        if (print_functional):
            print("--|" + self.topology_optimization_stage_str + "| ---> Vorticity Functional: (no weight)", self.fluid_functionals[2])    
        else:
            print("--|" + self.topology_optimization_stage_str + "| ---> Vorticity Functional")

    def _EvaluateOutletTransportScalarFunctional(self, print_functional=False):
        """
        This method computes the Outlet Concentration functional: int_{\Gamma_{out}}{1/2(c-c_target)^2}
        """
        if (self.first_iteration):
            self.SetTargetOutletTransportScalar()
        outlet_mp = self._FindAdjointSurfaceSourceProcess().model_part
        integral_value = 0.0
        integral_value_sq = 0.0
        for condition in outlet_mp.Conditions:
            geom = condition.GetGeometry()
            size = geom.DomainSize()  # Surface area in 3D
            nodes = condition.GetNodes()
            cond_transport_scalar = 0.0
            cond_transport_scalar_sq = 0.0
            for node in nodes:
                temp_target_diff = node.GetSolutionStepValue(KratosMultiphysics.TEMPERATURE) - self.target_outlet_transport_scalar
                cond_transport_scalar    += temp_target_diff
                cond_transport_scalar_sq += temp_target_diff**2
            cond_transport_scalar    /= len(nodes)
            cond_transport_scalar_sq /= len(nodes)
            integral_value    += cond_transport_scalar*size
            integral_value_sq += cond_transport_scalar_sq*size
        self.transport_functionals[0] = 0.5 * integral_value_sq
        self.avg_outlet_transport_scalar_diff = integral_value*size
        if (self.first_iteration):
            self.initial_transport_functionals_values[0] = self.transport_functionals[0] 
        if (print_functional):
            print("--|" + self.topology_optimization_stage_str + "| ---> Outlet Transport Scalar Functional (no weight):", self.transport_functionals[0])
        else:
            print("--|" + self.topology_optimization_stage_str + "| ---> Outlet Transport Scalar Functional")

    def _EvaluateFocusRegionTransportScalarFunctional(self, print_functional=False):
        """
        This method computes the Focus Region Concentration functional: int_{\Gamma_{out}}{c}
        """
        focus_mp = self._FindAdjointVolumeSourceProcess().model_part
        focus_nodes_list = [(node.Id-1) for node in focus_mp.Nodes]
        t_focus_sq = np.asarray(KratosMultiphysics.VariableUtils().GetSolutionStepValuesVector(focus_mp.Nodes, KratosMultiphysics.TEMPERATURE, 0))**2
        self.transport_functionals[1] = np.dot(self.nodal_domain_sizes[focus_nodes_list], t_focus_sq)
        if (self.first_iteration):
            self.initial_transport_functionals_values[1] = self.transport_functionals[1] 
        if (print_functional):
            print("--|" + self.topology_optimization_stage_str + "| ---> Focus Region Concentration Functional (no weight):", self.transport_functionals[1])
        else:
            print("--|" + self.topology_optimization_stage_str + "| ---> Focus Region Concentration Functional")

    def _EvaluateTransportScalarDiffusionFunctional(self, print_functional=False):
        """
        This method computes the Transport Scalar Diffusion functional: int_{\Omega}{D\\||grad(u)||^2}
        """
        mp = self._GetComputingModelPart()
        transport_scalar_gradient = np.asarray(KratosMultiphysics.VariableUtils().GetSolutionStepValuesVector(mp.Nodes, KratosMultiphysics.TEMPERATURE_GRADIENT, 0, self.dim)).reshape(self.n_nodes, self.dim)
        transport_scalar_gradient_norm_squared = (np.linalg.norm(transport_scalar_gradient, axis=1))**2
        integrand = self.conductivity * transport_scalar_gradient_norm_squared
        self.transport_functionals[2] = np.dot(integrand, self.nodal_domain_sizes)
        if (self.first_iteration):
            self.initial_transport_functionals_values[2] = self.transport_functionals[2]
        if (print_functional):
            print("--|" + self.topology_optimization_stage_str + "| ---> Transport Scalar Diffusion Functional (no weight):", self.transport_functionals[2])
        else:
            print("--|" + self.topology_optimization_stage_str + "| ---> Transport Scalar Diffusion Functional")

    def _EvaluateTransportScalarConvectionFunctional(self, print_functional=False):
        """
        This method computes the Transport Scalar Convection functional: int_{\Omega}{beta*T*dot(u,grad(T))}
        """
        mp = self._GetComputingModelPart()
        transport_scalar = np.asarray(KratosMultiphysics.VariableUtils().GetSolutionStepValuesVector(mp.Nodes, KratosMultiphysics.TEMPERATURE, 0))
        transport_scalar_gradient = np.asarray(KratosMultiphysics.VariableUtils().GetSolutionStepValuesVector(mp.Nodes, KratosMultiphysics.TEMPERATURE_GRADIENT, 0, self.dim)).reshape(self.n_nodes, self.dim)
        convection_velocity = np.asarray(KratosMultiphysics.VariableUtils().GetSolutionStepValuesVector(mp.Nodes, KratosMultiphysics.VELOCITY, 0, self.dim)).reshape(self.n_nodes, self.dim)
        integrand = self.convection_coefficient * transport_scalar * np.einsum('ij,ij->i', convection_velocity, transport_scalar_gradient)
        self.transport_functionals[3] = np.dot(integrand, self.nodal_domain_sizes)
        if (self.first_iteration):
            self.initial_transport_functionals_values[3] = self.transport_functionals[3]
        if (print_functional):
            print("--|" + self.topology_optimization_stage_str + "| ---> Transport Scalar Convection Functional (no weight):", self.transport_functionals[3])
        else:
            print("--|" + self.topology_optimization_stage_str + "| ---> Transport Scalar Convection Functional")

    def _EvaluateTransportScalarDecayFunctional(self, print_functional=False):
        """
        This method computes the Transport Scalar Decay functional: int_{\Omega}{kT^2}
        """
        mp = self._GetComputingModelPart()
        transport_scalar = np.asarray(KratosMultiphysics.VariableUtils().GetSolutionStepValuesVector(mp.Nodes, KratosMultiphysics.TEMPERATURE, 0))
        integrand = self.decay*(transport_scalar**2)
        self.transport_functionals[4] = np.dot(integrand, self.nodal_domain_sizes)
        if (self.first_iteration):
            self.initial_transport_functionals_values[4] = self.transport_functionals[4]
        if (print_functional):
            print("--|" + self.topology_optimization_stage_str + "| ---> Transport Scalar Decay Functional (no weight):", self.transport_functionals[4])
        else:
            print("--|" + self.topology_optimization_stage_str + "| ---> Transport Scalar Decay Functional")

    def _EvaluateTransportScalarSourceFunctional(self, print_functional=False):
        """
        This method computes the Transport Scalar Source functional: int_{\Omega}{-Q*T}
        """
        mp = self._GetComputingModelPart()
        transport_scalar = np.asarray(KratosMultiphysics.VariableUtils().GetSolutionStepValuesVector(mp.Nodes, KratosMultiphysics.TEMPERATURE, 0))
        transport_scalar_volume_flux = np.asarray(KratosMultiphysics.VariableUtils().GetSolutionStepValuesVector(mp.Nodes, KratosMultiphysics.HEAT_FLUX, 0))
        integrand = -transport_scalar_volume_flux*transport_scalar
        self.transport_functionals[5] = np.dot(integrand, self.nodal_domain_sizes)
        if (self.first_iteration):
            self.initial_transport_functionals_values[5] = self.transport_functionals[5]
        if (print_functional):
            print("--|" + self.topology_optimization_stage_str + "| ---> Transport Scalar Source Functional (no weight):", self.transport_functionals[5])
        else:
            print("--|" + self.topology_optimization_stage_str + "| ---> Transport Scalar Source Functional")

    def SetTargetOutletTransportScalar(self, target_transport_scalar=0.0):
        self.target_outlet_transport_scalar = target_transport_scalar

    def _ComputeFunctionalDerivativesFunctionalContribution(self):
        return self._ComputeFunctionalDerivativesFluidFunctionalContribution() + self._ComputeFunctionalDerivativesTransportFunctionalContribution()
    
    def _ComputeFunctionalDerivativesPhysicsContribution(self):
        return self._ComputeFunctionalDerivativesFluidPhysicsContribution() + self._ComputeFunctionalDerivativesTransportPhysicsContribution()

    def _GetFluidSolver(self):
        if (not self.IsAdjointStage()):
            return self._GetPhysicsSolver()._GetFluidSolver()
        else:
            return self._GetAdjointSolver()._GetFluidSolver()
    
    def _GetTransportSolver(self):
        if (not self.IsAdjointStage()):
            return self._GetPhysicsSolver()._GetTransportSolver()
        else:
            return self._GetAdjointSolver()._GetTransportSolver()
        
    def _GetFluidModelPart(self):
        return self._GetFluidSolver.GetMainModelPart()

    def _GetTransportModelPart(self):
        return self._GetTransportSolver.GetMainModelPart()

    def _CheckMaterialProperties(self):
        print("--|CHECK| Check Physics Properties")
        self._GetSolver()._CheckMaterialProperties()

    def _InitializePhysicsParameters(self):
        print("--|" + self.topology_optimization_stage_str + "| INITIALIZE PHYSICS PARAMETERS")
        self._InitializeResistance()
        self._InitializeConductivity()
        self._InitializeDecay()
        self._InitializeConvectionCoefficient()

    def ResetPhysicsParameters(self):
        self._ResetResistance()
        self._ResetConductivity()
        self._ResetDecay()
        self._ResetConvectionCoefficient()

    def UpdatePhysicsParametersVariables(self):
        self._UpdateResistanceVariable()
        self._UpdateConductivityVariable()
        self._UpdateDecayVariable()
        self._UpdateConvectionCoefficientVariable()

    def _UpdatePhysicsParameters(self):
        self._UpdateResistance()
        self._UpdateConductivity()
        self._UpdateDecay()
        self._UpdateConvectionCoefficient()
        




    








        