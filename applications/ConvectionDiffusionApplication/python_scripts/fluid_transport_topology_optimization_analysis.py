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

class FluidTransportTopologyOptimizationAnalysisTime(TransportTopologyOptimizationAnalysis):
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
    
    def _InitializeTopologyOptimizationProblem(self):
        super()._InitializeTopologyOptimizationProblem()        
        
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

    def _PrintPhysicsFunctionals(self):
        self._PrintFluidFunctionals()
        self._PrintTransportFunctionals()
        
    def _InitializeFunctionalWeights(self):
        "This method assumes that the '_ImportFunctionalWeights()' has already been called"
        if (self.functional_weights_imported):
            # get number of functionals
            self.n_fluid_functionals = len(self.normalized_fluid_functional_weights)
            self.n_transport_functionals = len(self.normalized_transport_functional_weights)
            self.n_functionals = self.n_fluid_functionals + self.n_transport_functionals
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
            self._SetNormalizationFunctionals()
            self.functional_weights = self._RescaleFunctionalWeightsByNormalizationValues()
        else:
            info_msg = "Calling '_InitializeFunctionalWeights' method before '_ImportFunctionalWeights()'"
            raise RuntimeError("FluidTopologyOptimizationAnalysis: " + info_msg)

    def _ImportFunctionalWeights(self):
        self._ImportCouplingFunctionalWeights()
        super()._ImportFunctionalWeights()

    def ImportPhysicsFunctionalWeights(self):
        self._ImportFluidFunctionalWeights()
        self._ImportTransportFunctionalWeights()

    def _ImportCouplingFunctionalWeights(self):
        coupling_weights = [0.0, 0.0]
        functional_weights_parameters = self.optimization_settings["optimization_problem_settings"]["functional_weights"]["coupling"]
        coupling_weights[0] = functional_weights_parameters["fluid"].GetDouble()
        coupling_weights[1] = functional_weights_parameters["transport"].GetDouble()
        self.normalized_coupling_functional_weights = self._NormalizeFunctionalWeights(np.asarray(coupling_weights))

    def _SetInitialFunctionals(self):
        super()._SetInitialFunctionals()
        self._SetInitialCouplingFunctionals()

    def _SetInitialPhysicsFunctionals(self):
        self._SetInitialFluidFunctionals()
        self._SetInitialTransportFunctionals()

    def _SetInitialCouplingFunctionals(self):
        self.initial_coupling_functionals = np.asarray([self.initial_fluid_functional, self.initial_transport_functional])
        self.initial_coupling_functionals_abs_value = np.abs(self.initial_coupling_functionals)

    def _SetNormalizationPhysicsFunctionals(self):
        self._SetNormalizationFluidFunctionals()
        self._SetNormalizationTransportFunctionals()
    
    def _RescaleFunctionalWeightsByNormalizationValues(self):
        fluid_functional_weights = self._RescaleFluidFunctionalWeightsByNormalizationValues()
        transport_functional_weights = self._RescaleTransportFunctionalWeightsByNormalizationValues()
        return np.concatenate((fluid_functional_weights, transport_functional_weights))
    
    def _RescaleFluidFunctionalWeightsByNormalizationValues(self):
        if (np.abs(self.normalized_coupling_functional_weights[0]) < 1e-10):
            fluid_functional_weights = np.zeros(self.normalized_fluid_functional_weights.size)
        else:
            fluid_functional_weights = self.normalized_coupling_functional_weights[0] * super()._RescaleFluidFunctionalWeightsByNormalizationValues()
        return fluid_functional_weights
    
    def _RescaleTransportFunctionalWeightsByNormalizationValues(self):
        if (np.abs(self.normalized_coupling_functional_weights[1]) < 1e-10):
            transport_functional_weights = np.zeros(self.normalized_transport_functional_weights.size)
        else:
            transport_functional_weights = self.normalized_coupling_functional_weights[1] * super()._RescaleTransportFunctionalWeightsByNormalizationValues()
        return transport_functional_weights

    def _PrintFunctionalWeightsPhysicsInfo(self):
        self.MpiPrint("--|" + self.topology_optimization_stage_str + "| INITIAL FLUID FUNCTIONAL: " + str(self.initial_fluid_functional))
        self.MpiPrint("--|" + self.topology_optimization_stage_str + "| INITIAL TRANSPORT FUNCTIONAL: " + str(self.initial_transport_functional))
        self.MpiPrint("--|" + self.topology_optimization_stage_str + "| NORMALIZED FLUID FUNCTIONAL WEIGHTS: " + str(self.normalized_fluid_functional_weights))
        self.MpiPrint("--|" + self.topology_optimization_stage_str + "| NORMALIZED TRANSPORT FUNCTIONAL WEIGHTS: " + str(self.normalized_transport_functional_weights))

    def _EvaluateOptimizationPhysicsRequiredGradients(self):
        self._EvaluateVelocityGradient()
        self._EvaluateTransportScalarGradient()

    def EvaluatePhysicsFunctionals(self, print_functional):
        self.EvaluateFluidFunctionals(print_functional)
        self.EvaluateTransportFunctionals(print_functional)

    def EvaluateTotalFunctional(self):
        self.functionals = np.concatenate((self.fluid_functionals, self.transport_functionals))
        self.weighted_functionals = self.functional_weights * self.functionals
        self.functional = np.sum(self.weighted_functionals)
        if (self.first_iteration):
            self.initial_functional = self.functional

    def _ComputeFunctionalDerivativesFunctionalContribution(self):
        fluid_contribution = self._ComputeFunctionalDerivativesFluidFunctionalContribution() 
        transport_contribution = self._ComputeFunctionalDerivativesTransportFunctionalContribution()
        return (fluid_contribution + transport_contribution)
    
    def _ComputeFunctionalDerivativesPhysicsContribution(self):
        fluid_contribution = self._ComputeFunctionalDerivativesFluidPhysicsContribution() 
        transport_contribution = self._ComputeFunctionalDerivativesTransportPhysicsContribution()
        return (fluid_contribution + transport_contribution)

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

    def InitializeAnalysisTimeSettings(self):
        self.InitializeCoupledAnalysisTimeStepsSettings()
        super().InitializeAnalysisTimeSettings()

    def InitializeCoupledAnalysisTimeStepsSettings(self):
        pass

    def _InitializePhysicsParameters(self):
        self._InitializeFluidPhysicsParameters()
        self._InitializeTransportPhysicsParameters()

    def _ResetPhysicsParameters(self):
        self._ResetFluidPhysicsParameters()
        self._ResetTransportPhysicsParameters()

    def _UpdatePhysicsParametersVariables(self):
        self._UpdateFluidPhysicsParametersVariables()
        self._UpdateTransportPhysicsParametersVariables()

    def _UpdatePhysicsParameters(self):
        self._UpdateFluidPhysicsParameters()
        self._UpdateTransportPhysicsParameters()

    def _UpdateRelevantPhysicsVariables(self):
        self._UpdateRelevantFluidPhysicsVariables()
        self._UpdateRelevantTransportPhysicsVariables()

    def _UpdateRelevantAdjointVariables(self):
        self._UpdateRelevantFluidAdjointVariables()
        self._UpdateRelevantTransportAdjointVariables()

    def _InitializeFunctionalsInTime(self):
        self._InitializeFluidFunctionalsInTime()
        self._InitializeTransportFunctionalsInTime()

    def _SynchronizePhysicsParametersVariables(self):
        self._SynchronizeFluidPhysicsParametersVariables()
        self._SynchronizeTransportPhysicsParametersVariables()

    def EvaluatePhysicsFunctionalInDeltaTime(self, time_step_id):
        self.EvaluateFluidFunctionalsInDeltaTime(time_step_id)
        self.EvaluateTransportFunctionalsInDeltaTime(time_step_id)

    def StoreTimeStepSolution(self):
        self.StoreTimeStepFluidSolution()
        self.StoreTimeStepTransportSolution()
    
    def GetDefaultOptimizationSettings(self):
        ##settings string in json format
        default_optimization_settings = KratosMultiphysics.Parameters("""
        {
            "optimization_problem_settings": {
                "functional_weights": {
                    "fluid_functionals": {
                    "normalization" : {
                        "type" : "initial",
                        "value": 0.0
                        },
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
                    "transport_functionals": {
                        "normalization" : {
                            "type" : "initial",
                            "value": 0.0
                            },
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
                        },
                        "decay_1st_order_transport" : {
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

    def ResetPhysicsTimeStepVariables(self):
        self._ResetFluidTimeStepVariables()
        self._ResetTransportTimeStepVariables()

    def _SetTimePhysicsSolutionStorage(self):
        self._SetTimeFluidSolutionStorage()
        self._SetTimeTransportSolutionStorage()

    def InitializePhysicsSolutionStep(self):
        self._InitializeFluidSolutionStep()
        self._InitializeTransportSolutionStep()

    def _InitializeTransportSolutionStep(self):
        # Since we are working in coupled fluid problems, for each time step the corrsponding velocity is pass to the transport physics in the solver SolveSolutionStep.
        # It is not necessary to pass it at each InitializeSolutionStep.
        # For this reason the following line has been commented w.r.t. the Transport Top Opt Analysis:
        # self._InitializeTransportSolutionStepConvectionVelocity()
        pass

    def InitializeAdjointPhysicsSolutionStep(self):
        self._InitializeAdjointFluidSolutionStep()
        self._InitializeAdjointTransportSolutionStep()

    def InitializePhysicsFunctionalsWeightsInTime(self):
        self._InitializeFluidFunctionalsWeightsInTime()
        self._InitializeTransportFunctionalsWeightsInTime()

    def _InitializeAdjointFluidSolutionStep(self):
        super()._InitializeAdjointFluidSolutionStep()
        # self._InitializeAdjointPhysicsSolutionStepFunctionalDerivativeTransportScalar()
        # Since we are working in coupled fluid problems, for each time step the corrsponding velocity is pass to the transport physics in the solver SolveSolutionStep.
        # It is not necessary to pass it at each InitializeSolutionStep.
        # For this reason the following line has been commented w.r.t. the Transport Top Opt Analysis:
        # self._InitializeAdjointPhysicsSolutionStepFunctionalDerivativeTransportScalarAdjoint()

    # def _InitializeAdjointPhysicsSolutionStepFunctionalDerivativeTransportScalarAdjoint(self):
    #     transport_scalar_adj = self.adjoint_transport_scalar_solutions[self.n_time_steps-self.time_step_counter,:].copy()
    #     self._GetAdjointSolver()._UpdateFunctionalDerivativeTransportScalarAdjointVariable(transport_scalar_adj)
    #     self._GetMainModelPart().GetCommunicator().SynchronizeNonHistoricalVariable(KratosMultiphysics.FUNCTIONAL_DERIVATIVE_TRANSPORT_SCALAR_ADJ)
        
        




    








        