from sys import argv

import numpy as np #import the numpy library
import scipy as sp #import the scipy library
from scipy.spatial import KDTree
from scipy.sparse import dok_matrix, lil_matrix
import mmapy as MMA # import the MMA subroutines python library: https://github.com/arjendeetman/GCMMA-MMA-Python
import time

import KratosMultiphysics as KratosMultiphysics
from KratosMultiphysics.kratos_utilities import IssueDeprecationWarning
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
from KratosMultiphysics.FluidDynamicsApplication import trilinos_fluid_topology_optimization_solver
from KratosMultiphysics.ConvectionDiffusionApplication import trilinos_transport_topology_optimization_solver
from KratosMultiphysics.ConvectionDiffusionApplication import trilinos_fluid_transport_topology_optimization_solver

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

class FluidTransportTopologyOptimizationAnalysis(TransportTopologyOptimizationAnalysis):
    def __init__(self,model,parameters):
        super().__init__(model,parameters) 

    def _SetTopologyOptimizationName(self):
        self.topology_optimization_name = "FLUID-TRANSPORT"

    def _CreateTopologyOptimizationSolvers(self):
        """
        This method creates the NS and ADJ_NS solvers
        """
        if self.IsMpiParallelism():
            self.physics_solver = trilinos_fluid_transport_topology_optimization_solver.CreateSolver(self.model, self.project_parameters)
            self.adjoint_solver = trilinos_fluid_transport_topology_optimization_solver.CreateSolver(self.model, self.project_parameters, isAdjointSolver=True)
        else:
            self.physics_solver = fluid_transport_topology_optimization_solver.CreateSolver(self.model, self.project_parameters)
            self.adjoint_solver = fluid_transport_topology_optimization_solver.CreateSolver(self.model, self.project_parameters, isAdjointSolver=True)

    def _CreateSolver(self, isAdjointSolver = False):
        """
        This method creates a solver
        isAdjointSolver == False --> physics_solver
        isAdjointSolver == True  --> adjoint_solver
        """
        if self.IsMpiParallelism():
            return trilinos_fluid_transport_topology_optimization_solver.CreateSolver(self.model, self.project_parameters, isAdjointSolver)
        else:
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
        problem_stage = 1: Navier Stokes + Transport solution
        problem_stage = 2: Adjoint Navier-Stokes + Adjoint Transport solution
        """
        super()._RunStageSolutionLoop()          
        
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
        self._PrintTotalFunctional()
        self._PrintFluidFunctionals()
        self._PrintTransportFunctionals()
        
    def _InitializeFunctionalWeights(self):
        fluid_functional_weights     = self._ImportFluidFunctionalWeights()
        transport_functional_weights = self._ImportTransportFunctionalWeights()
        coupling_functional_weights  = self._ImportCouplingFunctionalWeights()
        # normalize weights
        self.normalized_fluid_functional_weights = self._NormalizeFunctionalWeights(np.asarray(fluid_functional_weights))
        self.normalized_transport_functional_weights = self._NormalizeFunctionalWeights(np.asarray(transport_functional_weights))
        self.normalized_coupling_functional_weights = self._NormalizeFunctionalWeights(np.asarray(coupling_functional_weights))
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
            self.MpiPrint("[WARNING] Initial fluid functional is zero")
        if (abs(self.initial_transport_functional) < 1e-10):
            self.MpiPrint("[WARNING] Initial transport functional is zero")
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
    
    def _PrintFunctionalWeightsPhysicsInfo(self):
        self.MpiPrint("--|" + self.topology_optimization_stage_str + "| INITIAL FLUID FUNCTIONAL: " + str(self.initial_fluid_functional))
        self.MpiPrint("--|" + self.topology_optimization_stage_str + "| INITIAL TRANSPORT FUNCTIONAL: " + str(self.initial_transport_functional))
        self.MpiPrint("--|" + self.topology_optimization_stage_str + "| NORMALIZED FLUID FUNCTIONAL WEIGHTS: " + str(self.normalized_fluid_functional_weights))
        self.MpiPrint("--|" + self.topology_optimization_stage_str + "| NORMALIZED TRANSPORT FUNCTIONAL WEIGHTS: " + str(self.normalized_transport_functional_weights))

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

    def _EvaluateRequiredGradients(self):
        self._ComputeScalarVariableNodalGradient(KratosMultiphysics.DESIGN_PARAMETER, KratosMultiphysics.DESIGN_PARAMETER_GRADIENT)
        self._ComputeScalarVariableNodalGradient(KratosMultiphysics.VELOCITY_X, KratosMultiphysics.VELOCITY_X_GRADIENT)
        self._ComputeScalarVariableNodalGradient(KratosMultiphysics.VELOCITY_Y, KratosMultiphysics.VELOCITY_Y_GRADIENT)
        if (self.dim == 3):
            self._ComputeScalarVariableNodalGradient(KratosMultiphysics.VELOCITY_Z, KratosMultiphysics.VELOCITY_Z_GRADIENT)
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

    def MpiSynchronizeLocalFunctionalValues(self):
        local_values = self.functionals
        # Sum the values across all ranks
        total_values = self.data_communicator.SumAll(local_values)
        self.functionals = total_values
        if (self.MpiRunOnlyRank(0)):
            self.weighted_functionals  = self.functional_weights * self.functionals
            self.functional  = np.sum(self.weighted_functionals)

    def _EvaluateResistanceFunctional(self, print_functional=False):
        """
        This method computes the resistance functional: int_{\Omega}{\\alpha||u||^2}
        """
        mp = self._GetComputingModelPart()
        velocity = np.asarray(KratosMultiphysics.VariableUtils().GetSolutionStepValuesVector(self._GetLocalMeshNodes(mp), KratosMultiphysics.VELOCITY, 0, self.dim)).reshape(self.n_nodes, self.dim)
        nodal_velocity_norm = np.linalg.norm(velocity, axis=1)
        integrand = self.resistance * (nodal_velocity_norm**2) #component-wise multiplication
        self.fluid_functionals[0] = np.dot(self.nodal_domain_sizes, integrand)
        if _CheckIsDistributed():
            self.fluid_functionals[0] = self.MpiSumLocalValues(self.fluid_functionals[0])
        if (self.first_iteration):
            self.initial_fluid_functionals_values[0] = self.fluid_functionals[0] 
        if (print_functional):
            self.MpiPrint("--|" + self.topology_optimization_stage_str + "| ---> Resistance Functional (no weight): " + str(self.fluid_functionals[0]))
        else:
            self.MpiPrint("--|" + self.topology_optimization_stage_str + "| ---> Resistance Functional")

    def _EvaluateStrainRateFunctional(self, print_functional=False):
        """
        This method computes the Strain-Rate functional: int_{\Omega}{\\2*mu*||1/2*[grad(u)+grad(u)^T]||^2}
        """
        vel_gradient_on_nodes = self._AssembleVelocityGradientOnNodes()
        vel_symmetric_gradient = 1.0/2.0 * (vel_gradient_on_nodes+(np.transpose(vel_gradient_on_nodes, axes=(0,2,1))))
        vel_symmetric_gradient_norm_squared = (np.linalg.norm(vel_symmetric_gradient, ord='fro', axis=(1, 2)))**2
        mu = self._GetViscosity()
        self.fluid_functionals[1] = 2.0*mu* np.dot(vel_symmetric_gradient_norm_squared, self.nodal_domain_sizes)
        if _CheckIsDistributed():
            self.fluid_functionals[1] = self.MpiSumLocalValues(self.fluid_functionals[1])
        if (self.first_iteration):
            self.initial_fluid_functionals_values[1] = self.fluid_functionals[1]
        if (print_functional):
            self.MpiPrint("--|" + self.topology_optimization_stage_str + "| ---> Strain-Rate Functional (no weight): " + str(self.fluid_functionals[1]))
        else:
            self.MpiPrint("--|" + self.topology_optimization_stage_str + "| ---> Strain-Rate Functional")

    def _EvaluateVorticityFunctional(self, print_functional=False):
        """
        This method computes the Vorticity functional: int_{\Omega}{\\2*mu*||1/2*[grad(u)-grad(u)^T]||^2}
        """
        vel_gradient_on_nodes = self._AssembleVelocityGradientOnNodes()
        vel_antisymmetric_gradient = 0.5 * (vel_gradient_on_nodes-(np.transpose(vel_gradient_on_nodes, axes=(0,2,1))))
        vel_antisymmetric_gradient_norm_squared = (np.linalg.norm(vel_antisymmetric_gradient, ord='fro', axis=(1, 2)))**2
        mu = self._GetViscosity()
        self.fluid_functionals[2] = 2.0*mu* np.dot(vel_antisymmetric_gradient_norm_squared, self.nodal_domain_sizes)
        if _CheckIsDistributed():
            self.fluid_functionals[2] = self.MpiSumLocalValues(self.fluid_functionals[2])
        if (self.first_iteration):
            self.initial_fluid_functionals_values[2] = self.fluid_functionals[2]
        if (print_functional):
            self.MpiPrint("--|" + self.topology_optimization_stage_str + "| ---> Vorticity Functional: (no weight) " + str(self.fluid_functionals[2])) 
        else:
            self.MpiPrint("--|" + self.topology_optimization_stage_str + "| ---> Vorticity Functional")

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
        if _CheckIsDistributed():
            self.transport_functionals[0] = self.MpiSumLocalValues(self.transport_functionals[0])
        self.avg_outlet_transport_scalar_diff = integral_value*size
        if (self.first_iteration):
            self.initial_transport_functionals_values[0] = self.transport_functionals[0] 
        if (print_functional):
            self.MpiPrint("--|" + self.topology_optimization_stage_str + "| ---> Outlet Transport Scalar Functional (no weight): " + str(self.transport_functionals[0]))
        else:
            self.MpiPrint("--|" + self.topology_optimization_stage_str + "| ---> Outlet Transport Scalar Functional")

    def _EvaluateFocusRegionTransportScalarFunctional(self, print_functional=False):
        """
        This method computes the Focus Region Concentration functional: int_{\Gamma_{out}}{c}
        """
        if (self.first_iteration):
            self.SetTargetFocusRegionTransportScalar()
        t_target = self.target_focus_region_transport_scalar
        focus_mp = self._GetSubModelPart(self._GetTransportModelPart(), self.target_focus_region_transport_scalar_model_part_name)
        if focus_mp is None:
            focus_mp = self._GetOptimizationDomain()
        focus_nodes_list = [self.nodes_ids_global_to_local_partition_dictionary[node.Id] for node in self._GetLocalMeshNodes(focus_mp)]
        t_focus_sq = (np.asarray(KratosMultiphysics.VariableUtils().GetSolutionStepValuesVector(self._GetLocalMeshNodes(focus_mp), KratosMultiphysics.TEMPERATURE, 0))-t_target)**2
        self.transport_functionals[1] = np.dot(self.nodal_domain_sizes[focus_nodes_list], t_focus_sq)
        if _CheckIsDistributed():
            self.transport_functionals[1] = self.MpiSumLocalValues(self.transport_functionals[1])
        if (self.first_iteration):
            self.initial_transport_functionals_values[1] = self.transport_functionals[1] 
        if (print_functional):
            self.MpiPrint("--|" + self.topology_optimization_stage_str + "| ---> Focus Region Concentration Functional (no weight): " + str(self.transport_functionals[1]))
        else:
            self.MpiPrint("--|" + self.topology_optimization_stage_str + "| ---> Focus Region Concentration Functional")

    def _EvaluateTransportScalarDiffusionFunctional(self, print_functional=False):
        """
        This method computes the Transport Scalar Diffusion functional: int_{\Omega}{D\\||grad(u)||^2}
        """
        mp = self._GetComputingModelPart()
        transport_scalar_gradient = np.asarray(KratosMultiphysics.VariableUtils().GetSolutionStepValuesVector(self._GetLocalMeshNodes(mp), KratosMultiphysics.TEMPERATURE_GRADIENT, 0, self.dim)).reshape(self.n_nodes, self.dim)
        transport_scalar_gradient_norm_squared = (np.linalg.norm(transport_scalar_gradient, axis=1))**2
        integrand = self.conductivity * transport_scalar_gradient_norm_squared
        self.transport_functionals[2] = np.dot(integrand, self.nodal_domain_sizes)
        if _CheckIsDistributed():
            self.transport_functionals[2] = self.MpiSumLocalValues(self.transport_functionals[2])
        if (self.first_iteration):
            self.initial_transport_functionals_values[2] = self.transport_functionals[2]
        if (print_functional):
            self.MpiPrint("--|" + self.topology_optimization_stage_str + "| ---> Transport Scalar Diffusion Functional (no weight): " + str(self.transport_functionals[2]))
        else:
            self.MpiPrint("--|" + self.topology_optimization_stage_str + "| ---> Transport Scalar Diffusion Functional")

    def _EvaluateTransportScalarConvectionFunctional(self, print_functional=False):
        """
        This method computes the Transport Scalar Convection functional: int_{\Omega}{beta*T*dot(u,grad(T))}
        """
        mp = self._GetComputingModelPart()
        transport_scalar = np.asarray(KratosMultiphysics.VariableUtils().GetSolutionStepValuesVector(self._GetLocalMeshNodes(mp), KratosMultiphysics.TEMPERATURE, 0))
        transport_scalar_gradient = np.asarray(KratosMultiphysics.VariableUtils().GetSolutionStepValuesVector(self._GetLocalMeshNodes(mp), KratosMultiphysics.TEMPERATURE_GRADIENT, 0, self.dim)).reshape(self.n_nodes, self.dim)
        convection_velocity = np.asarray(KratosMultiphysics.VariableUtils().GetSolutionStepValuesVector(self._GetLocalMeshNodes(mp), KratosMultiphysics.VELOCITY, 0, self.dim)).reshape(self.n_nodes, self.dim)
        integrand = self.convection_coefficient * transport_scalar * np.einsum('ij,ij->i', convection_velocity, transport_scalar_gradient)
        self.transport_functionals[3] = np.dot(integrand, self.nodal_domain_sizes)
        if _CheckIsDistributed():
            self.transport_functionals[3] = self.MpiSumLocalValues(self.transport_functionals[3])
        if (self.first_iteration):
            self.initial_transport_functionals_values[3] = self.transport_functionals[3]
        if (print_functional):
            self.MpiPrint("--|" + self.topology_optimization_stage_str + "| ---> Transport Scalar Convection Functional (no weight): " + str(self.transport_functionals[3]))
        else:
            self.MpiPrint("--|" + self.topology_optimization_stage_str + "| ---> Transport Scalar Convection Functional")

    def _EvaluateTransportScalarDecayFunctional(self, print_functional=False):
        """
        This method computes the Transport Scalar Decay functional: int_{\Omega}{kT^2}
        """
        mp = self._GetComputingModelPart()
        transport_scalar = np.asarray(KratosMultiphysics.VariableUtils().GetSolutionStepValuesVector(self._GetLocalMeshNodes(mp), KratosMultiphysics.TEMPERATURE, 0))
        integrand = self.decay*(transport_scalar**2)
        self.transport_functionals[4] = np.dot(integrand, self.nodal_domain_sizes)
        if _CheckIsDistributed():
            self.transport_functionals[4] = self.MpiSumLocalValues(self.transport_functionals[4])
        if (self.first_iteration):
            self.initial_transport_functionals_values[4] = self.transport_functionals[4]
        if (print_functional):
            self.MpiPrint("--|" + self.topology_optimization_stage_str + "| ---> Transport Scalar Decay Functional (no weight): " + str(self.transport_functionals[4]))
        else:
            self.MpiPrint("--|" + self.topology_optimization_stage_str + "| ---> Transport Scalar Decay Functional")

    def _EvaluateTransportScalarSourceFunctional(self, print_functional=False):
        """
        This method computes the Transport Scalar Source functional: int_{\Omega}{-Q*T}
        """
        mp = self._GetComputingModelPart()
        transport_scalar = np.asarray(KratosMultiphysics.VariableUtils().GetSolutionStepValuesVector(self._GetLocalMeshNodes(mp), KratosMultiphysics.TEMPERATURE, 0))
        transport_scalar_volume_flux = np.asarray(KratosMultiphysics.VariableUtils().GetSolutionStepValuesVector(self._GetLocalMeshNodes(mp), KratosMultiphysics.HEAT_FLUX, 0))
        integrand = -transport_scalar_volume_flux*transport_scalar
        self.transport_functionals[5] = np.dot(integrand, self.nodal_domain_sizes)
        if _CheckIsDistributed():
            self.transport_functionals[5] = self.MpiSumLocalValues(self.transport_functionals[5])
        if (self.first_iteration):
            self.initial_transport_functionals_values[5] = self.transport_functionals[5]
        if (print_functional):
            self.MpiPrint("--|" + self.topology_optimization_stage_str + "| ---> Transport Scalar Source Functional (no weight): " + str(self.transport_functionals[5]))
        else:
            self.MpiPrint("--|" + self.topology_optimization_stage_str + "| ---> Transport Scalar Source Functional")

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
        return self._GetFluidSolver().GetMainModelPart()

    def _GetTransportModelPart(self):
        return self._GetTransportSolver().GetMainModelPart()

    def _CheckMaterialProperties(self, check = False):
        if (check):
            self.MpiPrint("--|CHECK| Check Physics Properties")
            self._GetSolver()._CheckMaterialProperties()

    def _InitializePhysicsParameters(self):
        self.MpiPrint("--|" + self.topology_optimization_stage_str + "| INITIALIZE PHYSICS PARAMETERS")
        self._InitializeResistance()
        self._InitializeConductivity()
        self._InitializeDecay()
        self._InitializeConvectionCoefficient()
        self._InitializeTransportSource()

    def ResetPhysicsParameters(self):
        self._ResetResistance()
        self._ResetConductivity()
        self._ResetDecay()
        self._ResetConvectionCoefficient()
        self._ResetTransportSource()

    def UpdatePhysicsParametersVariables(self):
        self._UpdateResistanceVariable()
        self._UpdateConductivityVariable()
        self._UpdateDecayVariable()
        self._UpdateConvectionCoefficientVariable()
        self._UpdateTransportSourceVariable()

    def _UpdatePhysicsParameters(self):
        self._UpdateResistance()
        self._UpdateConductivity()
        self._UpdateDecay()
        self._UpdateConvectionCoefficient()
        self._UpdateTransportSource()

    def _SynchronizePhysicsParametersVariables(self):
        # FLUID
        self._GetMainModelPart().GetCommunicator().SynchronizeNonHistoricalVariable(KratosCFD.RESISTANCE)
        self._GetMainModelPart().GetCommunicator().SynchronizeNonHistoricalVariable(KratosMultiphysics.CONDUCTIVITY)
        self._GetMainModelPart().GetCommunicator().SynchronizeNonHistoricalVariable(KratosCD.DECAY)
        self._GetMainModelPart().GetCommunicator().SynchronizeNonHistoricalVariable(KratosMultiphysics.CONVECTION_COEFFICIENT)
        self._GetMainModelPart().GetCommunicator().SynchronizeVariable(KratosMultiphysics.HEAT_FLUX)
    
    def GetDefaultOptimizationSettings(self):
        ##settings string in json format
        default_optimization_settings = KratosMultiphysics.Parameters("""
        {
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
        default_optimization_settings.AddMissingParameters(super().GetDefaultOptimizationSettings())
        return default_optimization_settings
        
        




    








        