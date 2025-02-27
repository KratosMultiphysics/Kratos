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
from KratosMultiphysics.analysis_stage import AnalysisStage
from KratosMultiphysics.FluidDynamicsApplication import fluid_topology_optimization_solver
from KratosMultiphysics.ConvectionDiffusionApplication import transport_topology_optimization_solver
from KratosMultiphysics.ConvectionDiffusionApplication import fluid_transport_topology_optimization_solver

class FluidTransportTopologyOptimizationAnalysis(FluidTopologyOptimizationAnalysis):
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
        self._DefineTransportProperties() must be done first to set the trnasport solver properties values before setting their relative Kratos VARIABLE
        """
        self._DefineTransportProperties()
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

    def _DefineTransportProperties(self, decay = 0.0, convection_coefficient = 1.0):
        self._DefineDecay(decay)
        self._DefineConvectionCoefficient(convection_coefficient)  
    
    def _DefineDecay(self, decay):
        self.decay = decay

    def _DefineConvectionCoefficient(self, convection_coefficient):
        self.convection_coefficient = convection_coefficient

    def _DefineConvectiveVelocity(self, const_convective_vel, velocity=[0,0,0]):
        self._GetTransportSolver()._DefineConvectionVelocity(const_convective_vel, velocity)
        self._GetAdjointTransportSolver()._DefineConvectionVelocity(const_convective_vel, velocity)
        
    def _GetPhysicsMainModelPartsList(self):
        return [self._GetPhysicsSolver()._GetFluidSolver().main_model_part, self._GetPhysicsSolver()._GetTransportSolver().main_model_part]

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
        if (process is not None):
            process.value += self.functional_weights[3]
        else:
            KratosMultiphysics.Logger.PrintError("'None' Adjoint Surface Source Process", "Trying to '_SetTransportSurfaceSourceFromFunctional' but the process doesn't exists.")

    def _FindAdjointSurfaceSourceProcess(self):
        # Get the list of processes
        process_list = self._GetListOfProcesses()
        # Find a process by class name
        for process in process_list:
            if isinstance(process, KratosMultiphysics.assign_scalar_variable_process.AssignScalarVariableProcess):
                if (process.variable == KratosCD.FACE_HEAT_FLUX_ADJ):
                    return process
        return None
    
    def _FindAdjointVolumeSourceProcess(self):
        # Get the list of processes
        process_list = self._GetListOfProcesses()
        # Find a process by class name
        for process in process_list:
            if isinstance(process, KratosMultiphysics.assign_scalar_variable_process.AssignScalarVariableProcess):
                if (process.variable == KratosCD.HEAT_FLUX_ADJ):
                    return process
        return None

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
            print("--|" + self.topology_optimization_stage_str + "| ---> Transport Scalar Transfer Functional (" + str(self.functional_weights[5]) + "):", self.weighted_functionals[5])

    def _InitializeFunctionalWeights(self, fluid_weights = [1, 1, 0], transport_weights=[1,0], coupling_weights=[1,1]):
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

    def _SetFunctionalWeights(self, fluid_weights = [1, 1, 0], transport_weights=[1,0], coupling_weights=[1,1]):
        print("--|" + self.topology_optimization_stage_str + "| ---> Set Functional Weights")
        self._InitializeFunctionalWeights(fluid_weights, transport_weights, coupling_weights)
        self.EvaluateFunctionals(print_functional=False)
        self._SetInitialFunctionals()
        self.functional_weights = self._RescaleFunctionalWeightsByInitialValues()
        self._GetComputingModelPart().ProcessInfo.SetValue(KratosMultiphysics.FUNCTIONAL_WEIGHTS, self.functional_weights)
        self.EvaluateTotalFunctional()
        self._SetTransportSurfaceSourceFromFunctional()

    def _SetInitialFunctionals(self):
        self.initial_fluid_functional = np.dot(self.normalized_fluid_functional_weights, self.initial_fluid_functionals_values)
        self.initial_transport_functional = np.dot(self.normalized_transport_functional_weights, self.initial_transport_functionals_values)
        self.initial_coupling_functionals = np.asarray([self.initial_fluid_functional, self.initial_transport_functional])
        self.initial_coupling_functionals_abs_value = np.abs(self.initial_coupling_functionals)

    def _RescaleFunctionalWeightsByInitialValues(self):
        # fluid
        fluid_functional_weights = self.normalized_fluid_functional_weights
        if (np.sum(np.abs(self.normalized_fluid_functional_weights)) > 1e-10):
            fluid_functional_weights *= self.normalized_coupling_functional_weights[0] / abs(self.initial_fluid_functional)
        # transport
        transport_functional_weights = self.normalized_transport_functional_weights
        if (np.sum(np.abs(self.normalized_transport_functional_weights)) > 1e-10):
            transport_functional_weights *= self.normalized_coupling_functional_weights[1] / abs(self.initial_transport_functional)
        return np.concatenate((fluid_functional_weights, transport_functional_weights))

    def _EvaluateFunctional(self, print_functional=False):
        """
        This method is used to evaluate the functional value
        # Functionals Database
        # 0: resistance  : int_{\Omega}{alpha*||u||^2}
        # 1: strain-rate : int_{\Omega}{2*mu*||S||^2} , with S = 1/2*(grad(u)+grad(u)^T) strain-rate tensor
        # 2: vorticity   : int_{\Omega}{2*mu*||R||^2} = int_{\Omega}{mu*||curl(u)||^2} , curl(u) = vorticity vector, R = 1/2*(grad(u)-grad(u)^T) rotation-rate tensor
        # 3: outlet_concentration : int_{\Gamma_{out}}{c}
        # 4: region_concentration: int_{\Omega}{c^2}
        # 4: ...
        # 5: ...
        # 5 is just a number big enough to contain the acutal database of functionals
        """
        self._SetTopologyOptimizationStage(3)
        print("--|" + self.topology_optimization_stage_str + "| EVALUATE FUNCTIONAL VALUE")
        self.EvaluateFunctionals(print_functional)
        self.EvaluateTotalFunctional()

    def EvaluateFunctionals(self, print_functional):
        mp = self._GetComputingModelPart()
        # FLUID FUNCTIONALS
        velocity = np.asarray(KratosMultiphysics.VariableUtils().GetSolutionStepValuesVector(mp.Nodes, KratosMultiphysics.VELOCITY, 0, self.dim)).reshape(self.n_nodes, self.dim)
        if (abs(self.normalized_fluid_functional_weights[0]) > 1e-10):
            self._EvaluateResistanceFunctional(velocity, print_functional)
        if (abs(self.normalized_fluid_functional_weights[1]) > 1e-10):
            self._EvaluateStrainRateFunctional(velocity, print_functional)
        if (abs(self.normalized_fluid_functional_weights[2]) > 1e-10):
            self._EvaluateVorticityFunctional(velocity, print_functional)
        # TRANSPORT FUNCTIONALS
        transport_scalar = np.asarray(KratosMultiphysics.VariableUtils().GetSolutionStepValuesVector(mp.Nodes, KratosMultiphysics.TEMPERATURE, 0))
        if (abs(self.normalized_transport_functional_weights[0]) > 1e-10):
            self._EvaluateOutletTransportScalarFunctional(print_functional)
        if (abs(self.normalized_transport_functional_weights[1]) > 1e-10):
            self._EvaluateFocusRegionTransportScalarFunctional(print_functional)
        if (abs(self.normalized_transport_functional_weights[2]) > 1e-10):
            self._EvaluateTransportScalarTransferFunctional(print_functional)

    def EvaluateTotalFunctional(self):
        self.functionals = np.concatenate((self.fluid_functionals, self.transport_functionals))
        self.weighted_functionals = self.functional_weights * self.functionals
        self.functional = np.sum(self.weighted_functionals)
        if (self.first_iteration):
            self.initial_functional = self.functional

    def _EvaluateResistanceFunctional(self, velocity, print_functional=False):
        """
        This method computes the resistance functional: int_{\Omega}{\\alpha||u||^2}
        """
        nodal_velocity_norm = np.linalg.norm(velocity, axis=1)
        integrand = self.resistance * (nodal_velocity_norm**2) #component-wise multiplication
        self.fluid_functionals[0] = np.dot(self.nodal_domain_sizes, integrand)
        if (self.first_iteration):
            self.initial_fluid_functionals_values[0] = self.fluid_functionals[0] 
        if (print_functional):
            print("--|" + self.topology_optimization_stage_str + "| ---> Resistance Functional (no weight):", self.functionals[0])
        else:
            print("--|" + self.topology_optimization_stage_str + "| --> Resistance Functional")

    def _EvaluateStrainRateFunctional(self, velocity, print_functional=False):
        """
        This method computes the Strain-Rate functional: int_{\Omega}{\\2*mu*||1/2*[grad(u)+grad(u)^T]||^2}
        """
        vel = velocity[self.element_nodes_ids[:]]
        vel_gradient = np.matmul(np.transpose(vel, axes=(0,2,1)), self.shape_functions_derivatives)
        vel_symmetric_gradient = 1.0/2.0 * (vel_gradient+(np.transpose(vel_gradient, axes=(0,2,1))))
        vel_symmetric_gradient_norm_squared = (np.linalg.norm(vel_symmetric_gradient, ord='fro', axis=(1, 2)))**2
        mu = self._GetComputingModelPart().Elements[1].Properties.GetValue(KratosMultiphysics.DYNAMIC_VISCOSITY)
        self.fluid_functionals[1] = 2*mu* np.dot(vel_symmetric_gradient_norm_squared, self.elemental_domain_size)
        if (self.first_iteration):
            self.initial_fluid_functionals_values[1] = self.fluid_functionals[1]
        if (print_functional):
            print("--|" + self.topology_optimization_stage_str + "| ---> Strain-Rate Functional (no weight):", self.functionals[1])
        else:
            print("--|" + self.topology_optimization_stage_str + "| --> Strain-Rate Functional")

    def _EvaluateVorticityFunctional(self, velocity, print_functional=False):
        """
        This method computes the Vorticity functional: int_{\Omega}{\\2*mu*||1/2*[grad(u)-grad(u)^T]||^2}
        """
        vel = velocity[self.element_nodes_ids[:]]
        vel_gradient = np.matmul(np.transpose(vel, axes=(0,2,1)), self.shape_functions_derivatives)
        vel_antisymmetric_gradient = 1.0/2.0 * (vel_gradient-(np.transpose(vel_gradient, axes=(0,2,1))))
        vel_antisymmetric_gradient_norm_squared = (np.linalg.norm(vel_antisymmetric_gradient, ord='fro', axis=(1, 2)))**2
        mu = self._GetComputingModelPart().Elements[1].Properties.GetValue(KratosMultiphysics.DYNAMIC_VISCOSITY)
        self.fluid_functionals[2] = 2*mu* np.dot(vel_antisymmetric_gradient_norm_squared, self.elemental_domain_size)
        if (self.first_iteration):
            self.initial_fluid_functionals_values[2] = self.fluid_functionals[2]
        if (print_functional):
            print("--|" + self.topology_optimization_stage_str + "| ---> Vorticity Functional: (no weight)", self.functionals[2])    
        else:
            print("--|" + self.topology_optimization_stage_str + "| --> Vorticity Functional")

    def _EvaluateOutletTransportScalarFunctional(self, print_functional=False):
        """
        This method computes the Outlet Concentration functional: int_{\Gamma_{out}}{c}
        """
        outlet_mp = self._FindAdjointSurfaceSourceProcess().model_part
        integral_value = 0.0
        for condition in outlet_mp.Conditions:
            geom = condition.GetGeometry()
            size = geom.DomainSize()  # Surface area in 3D
            nodes = condition.GetNodes()
            cond_transport_scalar = 0.0
            for node in nodes:
                cond_transport_scalar += node.GetSolutionStepValue(KratosMultiphysics.TEMPERATURE)**2
            cond_transport_scalar /= len(nodes)
            integral_value += cond_transport_scalar*size
        self.transport_functionals[0] = integral_value
        if (self.first_iteration):
            self.initial_transport_functionals_values[0] = self.transport_functionals[0] 
        if (print_functional):
            print("--|" + self.topology_optimization_stage_str + "| ---> Outlet Concentration Functional (no weight):", self.functionals[3])
        else:
            print("--|" + self.topology_optimization_stage_str + "| --> Outlet Concentration Functional")

    def _EvaluateFocusRegionTransportScalarFunctional(self, print_functional=False):
        """
        This method computes the Focus Region Concentration functional: int_{\Gamma_{out}}{c}
        """
        mp = self._GetMainModelPart()
        t_focus_sq = np.asarray(KratosMultiphysics.VariableUtils().GetSolutionStepValuesVector(mp.Nodes, KratosCD.OPTIMIZATION_TEMPERATURE, 0))**2
        self.transport_functionals[1] = np.dot(self.nodal_domain_sizes, t_focus_sq)
        # integral_value = 0.0
        # for elem in focus_mp.Elements:
        #     geom = elem.GetGeometry()
        #     size = geom.DomainSize()  # Surface area in 3D
        #     nodes = elem.GetNodes()
        #     elem_transport_scalar = 0.0
        #     for node in nodes:
        #         elem_transport_scalar += node.GetSolutionStepValue(KratosMultiphysics.TEMPERATURE)**2
        #     elem_transport_scalar /= len(nodes)
        #     integral_value += elem_transport_scalar*size
        # self.transport_functionals[1] = integral_value
        if (self.first_iteration):
            self.initial_transport_functionals_values[1] = self.transport_functionals[1] 
        if (print_functional):
            print("--|" + self.topology_optimization_stage_str + "| ---> Focus Region Concentration Functional (no weight):", self.functionals[4])
        else:
            print("--|" + self.topology_optimization_stage_str + "| --> Focus Region Concentration Functional")

    def _EvaluateTransportScalarTransferFunctional(self, print_functional=False):
        """
        This method computes the Strain-Rate functional: int_{\Omega}{\\2*mu*||1/2*[grad(u)+grad(u)^T]||^2}
        """
        mp = self._GetTransportModelPart()
        transport_scalar = np.asarray(KratosMultiphysics.VariableUtils().GetSolutionStepValuesVector(mp.Nodes, KratosMultiphysics.TEMPERATURE, 0))
        t = transport_scalar[self.element_nodes_ids[:]]
        t_gradient = np.einsum('ij,ijk->ik', t, self.shape_functions_derivatives)
        t_gradient_norm_squared = (np.linalg.norm(t_gradient, axis=1))**2
        conductivity = self._GetTransportModelPart().Elements[1].Properties.GetValue(KratosMultiphysics.CONDUCTIVITY)
        self.transport_functionals[2] = conductivity*np.dot(t_gradient_norm_squared, self.elemental_domain_size)
        if (self.first_iteration):
            self.initial_transport_functionals_values[2] = self.transport_functionals[2]
        if (print_functional):
            print("--|" + self.topology_optimization_stage_str + "| ---> Transport Scalar Transfer Functional (no weight):", self.functionals[5])
        else:
            print("--|" + self.topology_optimization_stage_str + "| --> Transport Scalar Transfer Functional")

    
    def _EvaluateFunctionalDerivatives(self, first_it = False):
        """
        This method is used to evaluate the functional derivatives w.r.t the design parameter
        """
        self._SetTopologyOptimizationStage(3)
        print("--|" + self.topology_optimization_stage_str + "| EVALUATE FUNCTIONAL DERIVATIVES")
        mp = self._GetComputingModelPart()
        velocity = np.asarray(KratosMultiphysics.VariableUtils().GetSolutionStepValuesVector(mp.Nodes, KratosMultiphysics.VELOCITY, 0, self.dim)).reshape(self.n_nodes, self.dim)
        velocity_adjoint= np.asarray(KratosMultiphysics.VariableUtils().GetSolutionStepValuesVector(mp.Nodes, KratosMultiphysics.VELOCITY_ADJ, 0, self.dim)).reshape(self.n_nodes, self.dim)
        temp_functional_derivatives_wrt_design = self.resistance_derivative_wrt_design * np.sum(velocity * (self.functional_weights[0]*velocity + velocity_adjoint), axis=1) * self.nodal_domain_sizes 
        self.functional_derivatives_wrt_design = temp_functional_derivatives_wrt_design
        for node in mp.Nodes:
            node.SetValue(KratosMultiphysics.FUNCTIONAL_DERIVATIVE, self.functional_derivatives_wrt_design[node.Id-1])
        if (first_it):
            self.initial_functional_derivatives_wrt_design = self.functional_derivatives_wrt_design
        self.functional_derivatives_wrt_design = np.asarray([self.functional_derivatives_wrt_design]).T

    def _InitializeTopologyOptimizationStepPhysicsSolution(self):
        super()._InitializeTopologyOptimizationStepPhysicsSolution()
        self._UpdateConductivity()

    def _UpdateConductivity(self):
        pass
    
    def _UpdateRelevantAdjointVariables(self):
        self._UptateOptimizationTemperatureVariable()

    def _UptateOptimizationTemperatureVariable(self):
        focus_mp = self._FindAdjointVolumeSourceProcess().model_part
        for node in focus_mp.Nodes:
            node.SetSolutionStepValue(KratosCD.OPTIMIZATION_TEMPERATURE, node.GetSolutionStepValue(KratosMultiphysics.TEMPERATURE))

    def _GetFluidModelPart(self):
        return self._GetPhysicsSolver()._GetFluidSolver().GetMainModelPart()

    def _GetTransportModelPart(self):
        return self._GetPhysicsSolver()._GetTransportSolver().GetMainModelPart()
    
    def PreparePhysicsSolver(self):
        super().PreparePhysicsSolver()
        self._SetDecay()
        self._SetConvection()
    
    def PrepareAdjointSolver(self):
        super().PrepareAdjointSolver()
        self._SetDecay()
        self._SetConvection()

    def _SetDecay(self):
        mp = self._GetTransportModelPart()
        for el in mp.Elements:
            el.Properties.SetValue(KratosCD.DECAY, self.decay)

    def _SetConvection(self):
        self._SetConvectiveVelocity()
        self._SetConvectionCoefficient()

    def _SetConvectionCoefficient(self):
        mp = self._GetTransportModelPart()
        for el in mp.Elements:
            el.Properties.SetValue(KratosMultiphysics.CONVECTION_COEFFICIENT, self.convection_coefficient)

    def _SetConvectiveVelocity(self):
        if (self._GetTransportSolver().IsConstantVelocity()):
            self._GetTransportSolver()._SetConstantConvectiveVelocity()
        else:
            self._GetTransportSolver()._SetNonConstantConvectiveVelocity()

    def _GetTransportSolver(self):
        if (self.IsAdjointStage()):
            return self._GetPhysicsSolver()._GetTransportSolver()
        else:
            return self._GetAdjointSolver()._GetTransportSolver()







    








        