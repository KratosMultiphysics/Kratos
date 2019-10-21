from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7

import KratosMultiphysics as KM
import KratosMultiphysics.StructuralMechanicsApplication.structural_analysis_for_thermal_coupling as structural_analysis_for_thermal_coupling
import KratosMultiphysics.StructuralMechanicsApplication.convection_diffussion_analysis_for_thermal_coupling as convection_diffussion_analysis_for_thermal_coupling

#============================================================================================================================
class MainThermoMechanicalAnalysis:
#============================================================================================================================

    def __init__(self, Model, StructuralParameters, ThermalParameters):
        self.MechanicalSolution = structural_analysis_for_thermal_coupling(Model, StructuralParameters)
        self.ThermalSolution = convection_diffussion_analysis_for_thermal_coupling(Model, ThermalParameters)

#============================================================================================================================

    def Run(self):
        self.Initialize()
        self.RunMainTemporalLoop()
        self.Finalize()

#============================================================================================================================

    def Initialize(self):
        self.ThermalSolution.Initialize()
        self.MechanicalSolution.Initialize()

#============================================================================================================================

    def RunMainTemporalLoop(self):
        # Temporal loop
        while self.MechanicalSolution.KeepAdvancingSolutionLoop():
            self.InitializeSolutionStep()
            self.SolveSolutionStep()
            self.FinalizeSolutionStep()

#============================================================================================================================

    def InitializeSolutionStep(self):
        self.ThermalSolution.time = self.ThermalSolution._GetSolver().AdvanceInTime(self.ThermalSolution.time)
        self.ThermalSolution.InitializeSolutionStep()
        self.MechanicalSolution.time = self.MechanicalSolution._GetSolver().AdvanceInTime(self.MechanicalSolution.time)
        self.MechanicalSolution.InitializeSolutionStep()

#============================================================================================================================

    def SolveSolutionStep(self):
            print("==================================================")
            print("==== Solving The Convection-Diffusion Problem ====")
            print("==================================================")
            self.ThermalSolution._GetSolver().Predict()
            is_converged = self.ThermalSolution._GetSolver().SolveSolutionStep()
            self.ThermalSolution.__CheckIfSolveSolutionStepReturnsAValue(is_converged)

            print("========================================")
            print("==== Solving The Mechanical Problem ====")
            print("========================================")
            self.MechanicalSolution._GetSolver().Predict()
            is_converged = self.MechanicalSolution._GetSolver().SolveSolutionStep()
            self.MechanicalSolution.__CheckIfSolveSolutionStepReturnsAValue(is_converged)
#============================================================================================================================
#============================================================================================================================
#============================================================================================================================
#============================================================================================================================
