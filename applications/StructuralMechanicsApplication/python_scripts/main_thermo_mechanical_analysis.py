from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7

import KratosMultiphysics as KM
import KratosMultiphysics.StructuralMechanicsApplication as SMA
import KratosMultiphysics.ConvectionDiffusionApplication as CDA
import KratosMultiphysics.StructuralMechanicsApplication.structural_analysis_for_thermal_coupling as structural_analysis_for_thermal_coupling
import KratosMultiphysics.StructuralMechanicsApplication.convection_diffussion_analysis_for_thermal_coupling as convection_diffussion_analysis_for_thermal_coupling

def Wait():
    input("Press Something")

#============================================================================================================================
class MainThermoMechanicalAnalysis:
#============================================================================================================================

    def __init__(self, Model, StructuralParameters, ThermalParameters):
        self.MechanicalSolution = structural_analysis_for_thermal_coupling.StructuralMechanicsAnalysisForThermalCoupling(Model, StructuralParameters)
        self.ThermalSolution = convection_diffussion_analysis_for_thermal_coupling.ConvectionDiffusionAnalysisForThermalCoupling(Model, ThermalParameters)

#============================================================================================================================

    def Run(self):
        self.Initialize()
        self.RunMainTemporalLoop()
        self.Finalize()

#============================================================================================================================

    def Initialize(self):
        self.ThermalSolution.Initialize()
        self.MechanicalSolution.Initialize()
        self.ShareNodesProcess()

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
            # self.ThermalSolution.__CheckIfSolveSolutionStepReturnsAValue(is_converged)

            print("========================================")
            print("==== Solving The Mechanical Problem ====")
            print("========================================")
            self.MechanicalSolution._GetSolver().Predict()
            is_converged = self.MechanicalSolution._GetSolver().SolveSolutionStep()
            # self.MechanicalSolution.__CheckIfSolveSolutionStepReturnsAValue(is_converged)

#============================================================================================================================

    def FinalizeSolutionStep(self):
        self.ThermalSolution.FinalizeSolutionStep()
        self.ThermalSolution.OutputSolutionStep()
        self.MechanicalSolution.FinalizeSolutionStep()
        self.MechanicalSolution.OutputSolutionStep()

#============================================================================================================================

    def Finalize(self):
        self.ThermalSolution.Finalize()
        self.MechanicalSolution.Finalize()

#============================================================================================================================

    def ShareNodesProcess(self):
        thermal_modelpart    = self.ThermalSolution.model.GetModelPart(self.ThermalSolution.project_parameters["solver_settings"]["model_part_name"].GetString())
        mechanical_modelpart = self.MechanicalSolution.model.GetModelPart(self.MechanicalSolution.project_parameters["solver_settings"]["model_part_name"].GetString())

        for node in thermal_modelpart.Nodes:
            node_id = node.Id
            node = mechanical_modelpart.GetNode(node_id)
        # self.ThermalSolution._GetSolver().AddVariables()
        # mechanical_modelpart.AddNodalSolutionStepVariable(KM.TEMPERATURE)
