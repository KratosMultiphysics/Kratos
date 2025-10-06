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
from KratosMultiphysics.ConvectionDiffusionApplication import transport_topology_optimization_solver
from KratosMultiphysics.ConvectionDiffusionApplication import trilinos_transport_topology_optimization_solver

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

class TransportTopologyOptimizationAnalysis(FluidTopologyOptimizationAnalysis):
    def _SetTopologyOptimizationName(self):
        self.topology_optimization_name = "TRANSPORT"
    
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
        N.B.: must be called after the creation of the transport model part
        problem_stage = 1: Navier Stokes solution
        problem_stage = 2: Adjoint Navier-Stokes solution
        """
        self._PrepareTransportSettings()
        super()._RunStageSolutionLoop()

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
            raise RuntimeError("Calling '_PrepareTransportSettings' in the wrong topology optimization stage.")
        
    def _GetPhysicsMainModelPartsList(self):
        return [self._GetPhysicsSolver().main_model_part]

    def PrintAnalysisStageProgressInformation(self):
        KratosMultiphysics.Logger.PrintInfo(self._GetSimulationName(), "TOTAL STEP: ", self._GetComputingModelPart().ProcessInfo[KratosMultiphysics.STEP])
        if self.IsPhysicsStage(): # NS
            KratosMultiphysics.Logger.PrintInfo(self._GetSimulationName(), " T STEP: ", self._GetComputingModelPart().ProcessInfo[KratosCD.TRANSPORT_TOP_OPT_T_STEP])
        elif self.IsAdjointStage(): # ADJ
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
            raise RuntimeError("Not Found Adjoint Volume Source Process", "It should be using the 'FACE_HEAT_FLUX_ADJ' variable")
    
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
            raise RuntimeError("Not Found Adjoint Volume Source Process", "It should be using the 'HEAT_FLUX_ADJ' variable")

    def _PrintFunctionals(self):
        self._PrintTotalFunctional()
        self._PrintTransportFunctionals()
    
    def _PrintTransportFunctionals(self):
        if (abs(self.functional_weights[3]) > 1e-10):
            self.MpiPrint("--|" + self.topology_optimization_stage_str + "| ---> Outlet Transport Scalar Functional (" + str(self.functional_weights[3]) + "): " + str(self.weighted_functionals[3]), min_echo=0)
        if (abs(self.functional_weights[4]) > 1e-10):
            self.MpiPrint("--|" + self.topology_optimization_stage_str + "| ---> Focus Region Transport Scalar Functional (" + str(self.functional_weights[4]) + "): " + str(self.weighted_functionals[4]), min_echo=0)
        if (abs(self.functional_weights[5]) > 1e-10):
            self.MpiPrint("--|" + self.topology_optimization_stage_str + "| ---> Transport Scalar Diffusion Functional (" + str(self.functional_weights[5]) + "): " + str(self.weighted_functionals[5]), min_echo=0)
        if (abs(self.functional_weights[6]) > 1e-10):
            self.MpiPrint("--|" + self.topology_optimization_stage_str + "| ---> Transport Scalar Convection Functional (" + str(self.functional_weights[6]) + "): " + str(self.weighted_functionals[6]), min_echo=0)
        if (abs(self.functional_weights[7]) > 1e-10):
            self.MpiPrint("--|" + self.topology_optimization_stage_str + "| ---> Transport Scalar Decay Functional (" + str(self.functional_weights[7]) + "): " + str(self.weighted_functionals[7]), min_echo=0)
        if (abs(self.functional_weights[8]) > 1e-10):
            self.MpiPrint("--|" + self.topology_optimization_stage_str + "| ---> Transport Scalar Source Functional (" + str(self.functional_weights[8]) + "): " + str(self.weighted_functionals[8]), min_echo=0)
        if (abs(self.functional_weights[9]) > 1e-10):
            self.MpiPrint("--|" + self.topology_optimization_stage_str + "| ---> Transport Scalar 1st Order Decay Functional (" + str(self.functional_weights[9]) + "): " + str(self.weighted_functionals[9]), min_echo=0)

    def _InitializeFunctionalWeights(self):
        transport_functional_weights = self._ImportTransportFunctionalWeights()
        # normalize weights
        self.normalized_transport_functional_weights = self._NormalizeFunctionalWeights(np.asarray(transport_functional_weights))
        # get number of functionals
        self.n_transport_functionals = len(self.normalized_transport_functional_weights)
        self.n_functionals = 3 + self.n_transport_functionals # 3 = transport fucntionals
        # initialize initial functionals vector container
        self.initial_transport_functionals_values = np.zeros(self.n_transport_functionals)
        self.initial_functionals_values = np.zeros(self.n_functionals)
        # initialize functionals vector container
        self.transport_functionals = np.zeros(self.n_transport_functionals)
        self.functionals = np.zeros(self.n_functionals)
        self.EvaluateFunctionals(print_functional=False)
        self._SetInitialFunctionals()
        self._SetNormalizationFunctionals()
        self.functional_weights = self._RescaleFunctionalWeightsByNormalizationValues()

    def _ImportTransportFunctionalWeights(self):
        transport_weights = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
        functional_weights_parameters = self.optimization_parameters["optimization_settings"]["optimization_problem_settings"]["functional_weights"]["transport_functionals"]
        self.transport_functional_normalization_strategy = functional_weights_parameters["normalization"]["type"].GetString()
        if self.transport_functional_normalization_strategy == "custom":
            self.transport_functional_normalization_value = functional_weights_parameters["normalization"]["value"].GetDouble()
        else:
            self.transport_functional_normalization_strategy == "initial"
        transport_weights[0] = functional_weights_parameters["outlet_transport_scalar"]["weight"].GetDouble()
        if (abs(transport_weights[0] > 1e-10)):
            raise RuntimeError("OUTLET_TRANSPORT_SCALAR FUNCTIONAL NOT WORKIUNG AND ITS WEIGHT IS DIFFERENT FROM ZERO", "Running '_ImportTransportFunctionalWeights' with the wrong transport_weights[0].")
        transport_weights[1] = functional_weights_parameters["focus_region_transport_scalar"]["weight"].GetDouble()
        transport_weights[2] = functional_weights_parameters["conductivity_transport"]["weight"].GetDouble()
        transport_weights[3] = functional_weights_parameters["convection_transport"]["weight"].GetDouble()
        transport_weights[4] = functional_weights_parameters["decay_transport"]["weight"].GetDouble()
        transport_weights[5] = functional_weights_parameters["source_transport"]["weight"].GetDouble()
        transport_weights[6] = functional_weights_parameters["decay_1st_order_transport"]["weight"].GetDouble()
        return transport_weights
    
    def _SetInitialFunctionals(self):
        self.initial_transport_functional = np.dot(self.normalized_transport_functional_weights, self.initial_transport_functionals_values)
        if (abs(self.initial_transport_functional) < 1e-10):
            self.MpiPrint("[WARNING] Initial transport functional is zero")

    def _SetNormalizationFunctionals(self):
        if self.transport_functional_normalization_strategy == "initial":
            self.transport_functional_normalization_value = self.initial_transport_functional
        else: # custom value, already defined
            pass

    def _RescaleFunctionalWeightsByNormalizationValues(self):
        if (np.sum(np.abs(self.normalized_transport_functional_weights)) < 1e-10):
            transport_functional_weights = np.zeros(self.normalized_transport_functional_weights.size)
        else:
            transport_functional_weights  = self.normalized_transport_functional_weights.copy()
            if (abs(self.initial_transport_functional) > 1e-15):
                transport_functional_weights /= abs(self.transport_functional_normalization_value)
        return np.concatenate((np.zeros(self.n_functionals-self.n_transport_functionals), transport_functional_weights))
    
    def _PrintFunctionalWeightsPhysicsInfo(self):
        self.MpiPrint("--|" + self.topology_optimization_stage_str + "| INITIAL FUNCTIONAL: " + str(self.initial_transport_functional))
        self.MpiPrint("--|" + self.topology_optimization_stage_str + "| NORMALIZED FUNCTIONAL WEIGHTS: " + str(self.normalized_transport_functional_weights))

    def _EvaluateRequiredGradients(self):
        self._ComputeScalarVariableNodalGradient(KratosMultiphysics.DESIGN_PARAMETER, KratosMultiphysics.DESIGN_PARAMETER_GRADIENT)
        self._ComputeScalarVariableNodalGradient(KratosMultiphysics.TEMPERATURE, KratosMultiphysics.TEMPERATURE_GRADIENT)
        self._ComputeScalarVariableNodalGradient(KratosMultiphysics.TEMPERATURE_ADJ, KratosMultiphysics.TEMPERATURE_ADJ_GRADIENT)

    def EvaluateFunctionals(self, print_functional):
        self._EvaluateRequiredGradients()
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
        if (abs(self.normalized_transport_functional_weights[6]) > 1e-10):
            self._EvaluateTransportScalar1stOrderDecayFunctional(print_functional)

    def EvaluateTotalFunctional(self):
        self.functionals = np.concatenate((np.zeros(self.n_functionals-self.n_transport_functionals), self.transport_functionals, ))
        self.weighted_functionals = self.functional_weights * self.functionals
        self.functional = np.sum(self.weighted_functionals)
        if (self.first_iteration):
            self.initial_functional = self.functional

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

    def _EvaluateTransportScalar1stOrderDecayFunctional(self, print_functional=False):
        """
        This method computes the Transport Scalar Decay functional: int_{\Omega}{kT^2}
        """
        mp = self._GetComputingModelPart()
        transport_scalar = np.asarray(KratosMultiphysics.VariableUtils().GetSolutionStepValuesVector(self._GetLocalMeshNodes(mp), KratosMultiphysics.TEMPERATURE, 0))
        integrand = self.decay*transport_scalar
        self.MpiPrint("integrand:"+integrand)
        self.transport_functionals[6] = np.dot(integrand, self.nodal_domain_sizes)
        if _CheckIsDistributed():
            self.transport_functionals[6] = self.MpiSumLocalValues(self.transport_functionals[6])
        if (self.first_iteration):
            self.initial_transport_functionals_values[6] = self.transport_functionals[6]
        if (print_functional):
            self.MpiPrint("--|" + self.topology_optimization_stage_str + "| ---> Transport Scalar 1st Order Decay Functional (no weight): " + str(self.transport_functionals[6]))
        else:
            self.MpiPrint("--|" + self.topology_optimization_stage_str + "| ---> Transport Scalar 1st Order Decay Functional")

    def SetTargetOutletTransportScalar(self):
        functional_settings = self.optimization_settings["optimization_problem_settings"]["functional_weights"]["transport_functionals"]["outlet_transport_scalar"]
        self.target_outlet_transport_scalar_model_part_name = functional_settings["outlet_model_part"].GetString()
        self.target_outlet_transport_scalar = functional_settings["target_value"].GetDouble()

    def SetTargetFocusRegionTransportScalar(self):
        functional_settings = self.optimization_settings["optimization_problem_settings"]["functional_weights"]["transport_functionals"]["focus_region_transport_scalar"]
        self.target_focus_region_transport_scalar_model_part_name = functional_settings["focus_region_model_part"].GetString()
        self.target_focus_region_transport_scalar = functional_settings["target_value"].GetDouble()

    def _ComputeFunctionalDerivativesFunctionalContribution(self):
        return self._ComputeFunctionalDerivativesTransportFunctionalContribution()
    
    def _ComputeFunctionalDerivativesPhysicsContribution(self):
        return self._ComputeFunctionalDerivativesTransportPhysicsContribution()

    def _ComputeFunctionalDerivativesTransportFunctionalContribution(self):
        local_mesh_nodes = self._GetLocalMeshNodes()
        velocity = np.asarray(KratosMultiphysics.VariableUtils().GetSolutionStepValuesVector(local_mesh_nodes, KratosMultiphysics.VELOCITY, 0, self.dim)).reshape(self.n_nodes, self.dim)
        transport_scalar = np.asarray(KratosMultiphysics.VariableUtils().GetSolutionStepValuesVector(local_mesh_nodes, KratosMultiphysics.TEMPERATURE, 0))
        transport_scalar_gradient = np.asarray(KratosMultiphysics.VariableUtils().GetSolutionStepValuesVector(local_mesh_nodes, KratosMultiphysics.TEMPERATURE_GRADIENT, 0, self.dim)).reshape(self.n_nodes, self.dim)
        transport_diffusion_functional_derivatives_wrt_design       = self.functional_weights[5] * self.conductivity_derivative_wrt_design_base * np.einsum('ij,ij->i', transport_scalar_gradient, transport_scalar_gradient) * self.nodal_domain_sizes
        transport_convection_functional_derivatives_wrt_design      = self.functional_weights[6] * self.convection_coefficient_derivative_wrt_design_base * (transport_scalar*(np.einsum('ij,ij->i', velocity, transport_scalar_gradient))) * self.nodal_domain_sizes
        transport_decay_functional_derivatives_wrt_design           = self.functional_weights[7] * self.decay_derivative_wrt_design_base * (transport_scalar**2) * self.nodal_domain_sizes
        transport_source_functional_derivatives_wrt_design          = self.functional_weights[8] * self.transport_source_derivative_wrt_design_base * transport_scalar * self.nodal_domain_sizes
        transport_1st_order_decay_functional_derivatives_wrt_design = self.functional_weights[9] * self.decay_derivative_wrt_design_base * transport_scalar * self.nodal_domain_sizes
        return transport_diffusion_functional_derivatives_wrt_design+transport_convection_functional_derivatives_wrt_design+transport_decay_functional_derivatives_wrt_design + transport_source_functional_derivatives_wrt_design + transport_1st_order_decay_functional_derivatives_wrt_design
    
    def _ComputeFunctionalDerivativesTransportPhysicsContribution(self):
        local_mesh_nodes = self._GetLocalMeshNodes()
        velocity = np.asarray(KratosMultiphysics.VariableUtils().GetSolutionStepValuesVector(local_mesh_nodes, KratosMultiphysics.VELOCITY, 0, self.dim)).reshape(self.n_nodes, self.dim)
        transport_scalar = np.asarray(KratosMultiphysics.VariableUtils().GetSolutionStepValuesVector(local_mesh_nodes, KratosMultiphysics.TEMPERATURE, 0))
        transport_scalar_adjoint = np.asarray(KratosMultiphysics.VariableUtils().GetSolutionStepValuesVector(local_mesh_nodes, KratosMultiphysics.TEMPERATURE_ADJ, 0))
        transport_scalar_gradient = np.asarray(KratosMultiphysics.VariableUtils().GetSolutionStepValuesVector(local_mesh_nodes, KratosMultiphysics.TEMPERATURE_GRADIENT, 0, self.dim)).reshape(self.n_nodes, self.dim)
        transport_scalar_adj_gradient = np.asarray(KratosMultiphysics.VariableUtils().GetSolutionStepValuesVector(local_mesh_nodes, KratosMultiphysics.TEMPERATURE_ADJ_GRADIENT, 0, self.dim)).reshape(self.n_nodes, self.dim)
        transport_physics_functional_derivatives_wrt_design  = self.conductivity_derivative_wrt_design_base * self.nodal_domain_sizes * np.einsum('ij,ij->i', transport_scalar_gradient, transport_scalar_adj_gradient)
        transport_physics_functional_derivatives_wrt_design += self.decay_derivative_wrt_design_base * (transport_scalar*transport_scalar_adjoint) * self.nodal_domain_sizes
        transport_physics_functional_derivatives_wrt_design += self.convection_coefficient_derivative_wrt_design_base * np.einsum('ij,ij->i', velocity, transport_scalar_gradient) * transport_scalar_adjoint * self.nodal_domain_sizes
        transport_physics_functional_derivatives_wrt_design += self.transport_source_derivative_wrt_design_base * transport_scalar_adjoint * self.nodal_domain_sizes
        return transport_physics_functional_derivatives_wrt_design
    
    def _UpdateRelevantAdjointVariables(self):
        super()._UpdateRelevantAdjointVariables()
        if (abs(self.functional_weights[3]) > 1e-10):
            if (self.first_iteration):
                self.SetTargetOutletTransportScalar()
            self._SetTransportSurfaceSourceFromFunctional()
        if (abs(self.functional_weights[4]) > 1e-10):
            if (self.first_iteration):
                self.SetTargetFocusRegionTransportScalar()
            self._UpdateOptimizationTemperatureVariable()

    def _UpdateOptimizationTemperatureVariable(self):
        focus_mp = self._GetSubModelPart(self._GetTransportModelPart(), self.target_focus_region_transport_scalar_model_part_name)
        if focus_mp is None:
            focus_mp = self._GetOptimizationDomain()
        target_t = self.target_focus_region_transport_scalar
        for node in self._GetLocalMeshNodes(focus_mp):
            node.SetSolutionStepValue(KratosCD.OPTIMIZATION_TEMPERATURE, node.GetSolutionStepValue(KratosMultiphysics.TEMPERATURE)-target_t)
        
    def _CheckMaterialProperties(self, check = False):
        if (check):
            self.MpiPrint("--|CHECK| Check Physics Properties")
            self._GetSolver()._CheckMaterialProperties()

    def _InitializePhysicsParameters(self):
        self.MpiPrint("--|" + self.topology_optimization_stage_str + "| INITIALIZE PHYSICS PARAMETERS")
        self._InitializeConductivity()
        self._InitializeDecay()
        self._InitializeConvectionCoefficient()
        self._InitializeTransportSource()
        self._SetConvectiveVelocity()    
    
    def _SetConvectiveVelocity(self, constant_velocity=False, vel=[0,0,0]):
        self.is_constant_velocity = constant_velocity
        self.constant_velocity = vel
        self._GetSolver()._SetConvectiveVelocity(self.is_constant_velocity, self.constant_velocity)

    def _InitializeConductivity(self):
        self.MpiPrint("--|" + self.topology_optimization_stage_str + "| INITIALIZE CONDUCTIVITY")
        self.conductivity_parameters = self.physics_parameters_settings["conductivity"] 
        self._ResetConductivity()
        self._UpdateConductivityVariable()

    def _InitializeDecay(self):
        self.MpiPrint("--|" + self.topology_optimization_stage_str + "| INITIALIZE DECAY")
        self.decay_parameters = self.physics_parameters_settings["decay"] 
        self._ResetDecay()
        self._UpdateDecayVariable()

    def _InitializeConvectionCoefficient(self):
        self.MpiPrint("--|" + self.topology_optimization_stage_str + "| INITIALIZE CONVECTION COEFFICIENT")
        self.convection_coefficient_parameters = self.physics_parameters_settings["convection_coefficient"] 
        self._ResetConvectionCoefficient()
        self._UpdateConvectionCoefficientVariable()

    def _InitializeTransportSource(self):
        self.MpiPrint("--|" + self.topology_optimization_stage_str + "| INITIALIZE TRANSPORT SOURCE")
        self.transport_source_parameters = self.physics_parameters_settings["transport_source"] 
        self._ResetTransportSource()
        self._UpdateTransportSourceVariable()

    def ResetPhysicsParameters(self):
        self._ResetConductivity()
        self._ResetDecay()
        self._ResetConvectionCoefficient()
        self._ResetTransportSource()

    def _ResetConductivity(self):
        self.conductivity = np.zeros(self.n_nodes)
        self.conductivity_derivative_wrt_design = np.zeros(self.n_nodes)

    def _ResetDecay(self):
        self.decay = np.zeros(self.n_nodes)
        self.decay_derivative_wrt_design = np.zeros(self.n_nodes)

    def _ResetConvectionCoefficient(self):
        self.convection_coefficient = np.zeros(self.n_nodes)
        self.convection_coefficient_derivative_wrt_design = np.zeros(self.n_nodes)

    def _ResetTransportSource(self):
        self.transport_source = np.zeros(self.n_nodes)
        self.transport_source_derivative_wrt_design = np.zeros(self.n_nodes)

    def UpdatePhysicsParametersVariables(self):
        self._UpdateConductivityVariable()
        self._UpdateDecayVariable()
        self._UpdateConvectionCoefficientVariable()
        self._UpdateTransportSourceVariable()

    def _UpdateConductivityVariable(self):
        self._GetSolver()._UpdateConductivityVariable(self.conductivity)

    def _UpdateDecayVariable(self):
        self._GetSolver()._UpdateDecayVariable(self.decay)

    def _UpdateConvectionCoefficientVariable(self):
        self._GetSolver()._UpdateConvectionCoefficientVariable(self.convection_coefficient)

    def _UpdateTransportSourceVariable(self):
        self._GetSolver()._UpdateTransportSourceVariable(self.transport_source)

    def _UpdatePhysicsParameters(self):
        self._UpdateConductivity()
        self._UpdateDecay()
        self._UpdateConvectionCoefficient()
        self._UpdateTransportSource()

    def _UpdateConductivity(self):
        """
        This method handles the conductivity update.
        """
        self.MpiPrint("--|" + self.topology_optimization_stage_str + "| UPDATE CONDUCTIVITY")
        self.conductivity, self.conductivity_derivative_wrt_design_base = self._ComputeConductivity(self.design_parameter)
        self._UpdateConductivityDesignDerivative()
        self._UpdateConductivityVariable()

    def _UpdateDecay(self):
        """
        This method handles the decay update.
        """
        self.MpiPrint("--|" + self.topology_optimization_stage_str + "| UPDATE DECAY")
        self.decay, self.decay_derivative_wrt_design_base = self._ComputeDecay(self.design_parameter)
        self._UpdateDecayDesignDerivative()
        self._UpdateDecayVariable()

    def _UpdateConvectionCoefficient(self):
        """
        This method handles the convection coefficient update.
        """
        self.MpiPrint("--|" + self.topology_optimization_stage_str + "| UPDATE CONVECTION COEFFICIENT")
        self.convection_coefficient, self.convection_coefficient_derivative_wrt_design_base = self._ComputeConvectionCoefficient(self.design_parameter)
        self._UpdateConvectionCoefficientDesignDerivative()
        self._UpdateConvectionCoefficientVariable()

    def _UpdateTransportSource(self):
        """
        This method handles the convection coefficient update.
        """
        self.MpiPrint("--|" + self.topology_optimization_stage_str + "| UPDATE TRANSPORT SOURCE")
        self.transport_source, self.transport_source_derivative_wrt_design_base = self._ComputeTransportSource(self.design_parameter)
        self._UpdateTransportSourceDesignDerivative()
        self._UpdateTransportSourceVariable()
    
    def _ComputeConductivity(self, design_parameter):
        return self._ComputePhysicsParameter(self.conductivity_parameters, design_parameter)
    
    def _ComputeDecay(self, design_parameter):
        return self._ComputePhysicsParameter(self.decay_parameters, design_parameter)
    
    def _ComputeConvectionCoefficient(self, design_parameter):
        return self._ComputePhysicsParameter(self.convection_coefficient_parameters, design_parameter)
    
    def _ComputeTransportSource(self, design_parameter):
        return self._ComputePhysicsParameter(self.transport_source_parameters, design_parameter)

    def _UpdateConductivityDesignDerivative(self):
        mask = self._GetOptimizationDomainNodesMask()
        conductivity_derivative_wrt_design_projected = self.conductivity_derivative_wrt_design_base * self.design_parameter_projected_derivatives
        self.conductivity_derivative_wrt_design = conductivity_derivative_wrt_design_projected
        # self.conductivity_derivative_wrt_design[mask] = self._ApplyDiffusiveFilterDerivative(conductivity_derivative_wrt_design_projected[mask])

    def _UpdateDecayDesignDerivative(self):
        mask = self._GetOptimizationDomainNodesMask()
        decay_derivative_wrt_design_projected = self.decay_derivative_wrt_design_base * self.design_parameter_projected_derivatives
        self.decay_derivative_wrt_design = decay_derivative_wrt_design_projected
        # self.decay_derivative_wrt_design[mask] = self._ApplyDiffusiveFilterDerivative(decay_derivative_wrt_design_projected[mask])

    def _UpdateConvectionCoefficientDesignDerivative(self):
        mask = self._GetOptimizationDomainNodesMask()
        convection_coefficient_derivative_wrt_design_projected = self.convection_coefficient_derivative_wrt_design_base * self.design_parameter_projected_derivatives
        self.convection_coefficient_derivative_wrt_design = convection_coefficient_derivative_wrt_design_projected
        # self.convection_coefficient_derivative_wrt_design[mask] = self._ApplyDiffusiveFilterDerivative(convection_coefficient_derivative_wrt_design_projected[mask])

    def _UpdateTransportSourceDesignDerivative(self):
        mask = self._GetOptimizationDomainNodesMask()
        transport_source_derivative_wrt_design_projected = self.transport_source_derivative_wrt_design_base * self.design_parameter_projected_derivatives
        self.transport_source_derivative_wrt_design = transport_source_derivative_wrt_design_projected
        # self.transport_source_derivative_wrt_design[mask] = self._ApplyDiffusiveFilterDerivative(transport_source_derivative_wrt_design_projected[mask])

    def _GetTransportModelPart(self):
        return self._GetMainModelPart()

    def GetDefaultPhysicsParametersSettings(self):
        ##settings string in json format
        default_physics_parameters_settings = KratosMultiphysics.Parameters("""
        {
            "conductivity": [{
                "interpolation_method": "polynomial",
                "value_void"        : 1e-4,
                "value_full"        : 1e-4,
                "interpolation_slope": {
                    "initial_slope": 1.0,
                    "final_slope"  : 10.0,
                    "iterations"   : [2,50]
                },
                "domain": "all"
            }],
            "decay": [{
                "interpolation_method": "polynomial",
                "value_void"        : 0.0,
                "value_full"        : 0.0,
                "interpolation_slope": {
                    "initial_slope": 1.0,
                    "final_slope"  : 10.0,
                    "iterations"   : [2,50]
                },
                "domain": "all"
            }],
            "convection_coefficient": [{
                "interpolation_method": "polynomial",
                "value_void"        : 1.0,
                "value_full"        : 1.0,
                "interpolation_slope": {
                    "initial_slope": 1.0,
                    "final_slope"  : 10.0,
                    "iterations"   : [2,50]
                },
                "domain": "all"
            }],
            "transport_source": [{
                "interpolation_method": "polynomial",
                "value_void"        : 0.0,
                "value_full"        : 0.0,
                "interpolation_slope": {
                    "initial_slope": 1.0,
                    "final_slope"  : 10.0,
                    "iterations"   : [2,50]
                },
                "domain": "all"
            }]
        }""")
        default_physics_parameters_settings.AddMissingParameters(super().GetDefaultPhysicsParametersSettings())
        return default_physics_parameters_settings
    
    def GetDefaultOptimizationSettings(self):
        ##settings string in json format
        default_optimization_settings = KratosMultiphysics.Parameters("""
        {
            "optimization_problem_settings": {
                "functional_weights": {
                    "normalization" : {
                        "type" : "initial",
                        "value": 0.0
                        },
                    "transport_functionals": {
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
                        "conductivity_transport" : {
                            "weight": 1.0
                        },
                        "convection_transport" : {
                            "weight": 0.0
                        },
                        "decay_transport" : {
                            "weight": 0.0
                        },
                        "source_transport" : {
                            "weight": 0.0
                        },
                        "decay_1st_order_transport" : {
                            "weight": 0.0
                        }
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

    def _CreateTopologyOptimizationSolvers(self):
        """
        This method creates the NS and ADJ_NS solvers
        """
        if self.IsMpiParallelism():
            self.physics_solver = trilinos_transport_topology_optimization_solver.CreateSolver(self.model, self.project_parameters)
            self.adjoint_solver = trilinos_transport_topology_optimization_solver.CreateSolver(self.model, self.project_parameters, isAdjointSolver=True)
        else:
            self.physics_solver = transport_topology_optimization_solver.CreateSolver(self.model, self.project_parameters)
            self.adjoint_solver = transport_topology_optimization_solver.CreateSolver(self.model, self.project_parameters, isAdjointSolver=True)

    def _CreateSolver(self, isAdjointSolver = False):
        """
        This method creates a solver
        isAdjointSolver == False --> physics_solver
        isAdjointSolver == True  --> adjoint_solver
        """
        if self.IsMpiParallelism():
            return trilinos_transport_topology_optimization_solver.CreateSolver(self.model, self.project_parameters, isAdjointSolver)
        else:
            return transport_topology_optimization_solver.CreateSolver(self.model, self.project_parameters, isAdjointSolver)
    
    def _SynchronizePhysicsParametersVariables(self):
        self._GetMainModelPart().GetCommunicator().SynchronizeNonHistoricalVariable(KratosMultiphysics.CONDUCTIVITY)
        self._GetMainModelPart().GetCommunicator().SynchronizeNonHistoricalVariable(KratosCD.DECAY)
        self._GetMainModelPart().GetCommunicator().SynchronizeNonHistoricalVariable(KratosMultiphysics.CONVECTION_COEFFICIENT)
        self._GetMainModelPart().GetCommunicator().SynchronizeVariable(KratosMultiphysics.HEAT_FLUX)
    
    

    


