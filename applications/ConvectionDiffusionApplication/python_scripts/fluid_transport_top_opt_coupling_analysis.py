
# Importing Kratos
import KratosMultiphysics
from KratosMultiphysics.ConvectionDiffusionApplication import python_solvers_wrapper_convection_diffusion as solver_wrapper

# Importing the base class
from KratosMultiphysics.analysis_stage import AnalysisStage
import KratosMultiphysics.ConvectionDiffusionApplication as KratosCD
import KratosMultiphysics.FluidDynamicsApplication as KratosCFD
from KratosMultiphysics.ConvectionDiffusionApplication.convection_diffusion_analysis import ConvectionDiffusionAnalysis
from KratosMultiphysics.ConvectionDiffusionApplication.transport_topology_optimization_analysis_test import TransportTopologyOptimizationAnalysisTest
from KratosMultiphysics.ConvectionDiffusionApplication import fluid_transport_topology_optimization_solver

# Other imports
import sys

class FluidTransportTopOptCouplingAnalysis(ConvectionDiffusionAnalysis):
    def __init__(self, model, project_parameters):
        super().__init__(model, project_parameters)
        self.topology_optimization_stage_str = "TEST"
    
    #### Internal functions ####
    def _CreateSolver(self):
        """ Create the Solver (and create and import the ModelPart if it is not alread in the model) """
        ## Solver construction
        return fluid_transport_topology_optimization_solver.CreateSolver(self.model, self.project_parameters)
    
    def InitializeSolutionStep(self):
        """This function performs all the required operations that should be executed
        (for each step) BEFORE solving the solution step.
        """
        self._GetSolver()._SetTopologyOptimizationStage(1)
        super().InitializeSolutionStep()

    def _GetFluidSolver(self):
        return self._GetSolver()._GetFluidSolver()

    def _GetTransportSolver(self):
        return self._GetSolver()._GetTransportSolver()
    
    def RunSolutionLoop(self):
        self.first_iteration = True
        self._InitializePhysicsParameters()
        self._InitializeTopologyOptimizationStepPhysicsSolution()
        self._CheckMaterialProperties()
        while self.KeepAdvancingSolutionLoop():
            self.time = self._AdvanceTime()
            self.InitializeSolutionStep()
            self.UpdatePhysicsParametersVariables()
            self._GetSolver().Predict()
            is_converged = self._GetSolver().SolveSolutionStep()
            # self.__CheckIfSolveSolutionStepReturnsAValue(is_converged)
            self.FinalizeSolutionStep()
            self.OutputSolutionStep()

    def _CheckMaterialProperties(self):
        print("--|CHECK| Check Physics Properties")
        self._GetSolver()._CheckMaterialProperties()

    def _InitializePhysicsParameters(self, resistance_parameters=[0.0, 1e4, 1.0], conductivity_parameters=[1e-4, 1e-2, 1.0], decay_parameters=[0.0, 0.0, 1.0], convection_coefficient_parameters=[0.0, 0.0, 1.0]):
        print("--|" + self.topology_optimization_stage_str + "| INITIALIZE PHYSICS PARAMETERS")
        self._InitializeResistance(resistance_parameters)
        self._InitializeConductivity(conductivity_parameters)
        self._InitializeDecay(decay_parameters)
        self._InitializeConvectionCoefficient(convection_coefficient_parameters)

    def _InitializeResistance(self, resistance_parameters_values):
        print("--|" + self.topology_optimization_stage_str + "| INITIALIZE RESISTANCE")
        self.resistance_parameters = resistance_parameters_values
        self._ResetResistance()
        self._UpdateResistanceVariable()

    def _InitializeConductivity(self, conductivity_parameters_values):
        print("--|" + self.topology_optimization_stage_str + "| INITIALIZE CONDUCTIVITY")
        self.conductivity_parameters = conductivity_parameters_values
        self._ResetConductivity()
        self._UpdateConductivityVariable()

    def _InitializeDecay(self, decay_parameters_values):
        print("--|" + self.topology_optimization_stage_str + "| INITIALIZE DECAY")
        self.decay_parameters = decay_parameters_values
        self._ResetDecay()
        self._UpdateDecayVariable()

    def _InitializeConvectionCoefficient(self, convection_coefficient_parameters_values):
        print("--|" + self.topology_optimization_stage_str + "| INITIALIZE CONVECTION COEFFICIENT")
        self.convection_coefficient_parameters = convection_coefficient_parameters_values
        self._ResetConvectionCoefficient()
        self._UpdateConvectionCoefficientVariable()

    def ResetPhysicsParameters(self):
        self._ResetResistance()
        self._ResetConductivity()
        self._ResetDecay()
        self._ResetConvectionCoefficient()

    def _ResetResistance(self):
        self.resistance = 0.0

    def _ResetConductivity(self):
        self.conductivity = 0.0

    def _ResetDecay(self):
        self.decay = 0.0

    def _ResetConvectionCoefficient(self):
        self.convection_coefficient = 0.0

    def UpdatePhysicsParametersVariables(self):
        self._UpdateResistanceVariable()
        self._UpdateConductivityVariable()
        self._UpdateDecayVariable()
        self._UpdateConvectionCoefficientVariable()

    def _UpdateResistanceVariable(self):
        self._GetSolver()._UpdateResistanceVariable(self.resistance)

    def _UpdateConductivityVariable(self):
        self._GetSolver()._UpdateConductivityVariable(self.conductivity)

    def _UpdateDecayVariable(self):
        self._GetSolver()._UpdateDecayVariable(self.decay)

    def _UpdateConvectionCoefficientVariable(self):
        self._GetSolver()._UpdateConvectionCoefficientVariable(self.convection_coefficient)

    def _UpdatePhysicsParameters(self):
        self._UpdateResistance()
        self._UpdateConductivity()
        self._UpdateDecay()
        self._UpdateConvectionCoefficient()
        
    def _UpdateResistance(self):
        """
        This method handles the resistance update.
        """
        print("--|" + self.topology_optimization_stage_str + "| UPDATE RESISTANCE")
        self.resistance = self._ComputeResistance()
        self._UpdateResistanceVariable()

    def _UpdateConductivity(self):
        """
        This method handles the conductivity update.
        """
        print("--|" + self.topology_optimization_stage_str + "| UPDATE CONDUCTIVITY")
        self.conductivity = self._ComputeConductivity()
        self._UpdateConductivityVariable()

    def _UpdateDecay(self):
        """
        This method handles the decay update.
        """
        print("--|" + self.topology_optimization_stage_str + "| UPDATE DECAY")
        self.decay = self._ComputeDecay()
        self._UpdateDecayVariable()

    def _UpdateConvectionCoefficient(self):
        """
        This method handles the condcutivity update.
        """
        print("--|" + self.topology_optimization_stage_str + "| UPDATE CONVECTION COEFFICIENT")
        self.convection_coefficient = self._ComputeConvectionCoefficient()
        self._UpdateConvectionCoefficientVariable()
        
    def _ComputeResistance(self):
        resistance = self.resistance_parameters[0]
        return resistance
    
    def _ComputeConductivity(self):
        conductivity = self.conductivity_parameters[0]
        return conductivity
    
    def _ComputeDecay(self):
        decay = self.decay_parameters[0]
        return decay
    
    def _ComputeConvectionCoefficient(self):
        convection_coeff = self.convection_coefficient_parameters[0]
        return convection_coeff

    def _InitializeTopologyOptimizationStepPhysicsSolution(self):
        """
        This method Re-Initializes the physical values in order to be consistent with the solution of the current optimization step
        """
        print("\n--|INITIALIZE OPTIMIZATION STEP PHYSICS SOLUTION|")
        self._GetSolver()._SetTopologyOptimizationStage(0)
        if (not self.first_iteration):
                self._ReInitializePhysics()
        self._UpdatePhysicsParameters()

    def _ReInitializePhysics(self):
        """
        This method Re-Initializes the physical values in order to be consistent with the solution of the current optimization step
        """
        print("--|" + self.topology_optimization_stage_str + "| RE-INITIALIZE PHYSICS")
        self._GetComputingModelPart().ProcessInfo.SetValue(KratosMultiphysics.TIME, 0.0)


    
    

    
