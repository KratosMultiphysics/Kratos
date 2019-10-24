from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
# importing the Kratos Library
import KratosMultiphysics
import KratosMultiphysics.CompressiblePotentialFlowApplication as KCPFApp
from KratosMultiphysics.CompressiblePotentialFlowApplication.potential_flow_adjoint_solver import PotentialFlowAdjointSolver
from KratosMultiphysics.CompressiblePotentialFlowApplication.potential_flow_solver import PotentialFlowFormulation

def CreateSolver(model, custom_settings):
    return PotentialFlowStochasticAdjointSolver(model, custom_settings)

class PotentialFlowStochasticAdjointSolver(PotentialFlowAdjointSolver):
    def __init__(self, model, custom_settings):

        # Construct the parent solver.
        super(PotentialFlowStochasticAdjointSolver, self).__init__(model, custom_settings)

        KratosMultiphysics.Logger.PrintInfo("::[PotentialFlowStochasticAdjointSolver]:: ", "Construction finished")

    def SolveSolutionStep(self):


        cvar_beta = self.response_function_settings["cvar_beta"].GetDouble()
        cvar_t = self.response_function_settings["cvar_t"].GetDouble()

        objective_function = self.response_function.CalculateValue(self.main_model_part)

        ##USING ABS
        if abs(objective_function) >= abs(cvar_t):
            multiplier = 1.0 / (1.0-cvar_beta)
            super(PotentialFlowStochasticAdjointSolver, self).SolveSolutionStep()
        else:
            multiplier = 0.0

        for node in self.main_model_part.Nodes:
            current_value = node.GetSolutionStepValue(KCPFApp.ADJOINT_VELOCITY_POTENTIAL)
            current_value_aux = node.GetSolutionStepValue(KCPFApp.ADJOINT_AUXILIARY_VELOCITY_POTENTIAL)
            node.SetSolutionStepValue(KCPFApp.ADJOINT_VELOCITY_POTENTIAL, multiplier*current_value)
            node.SetSolutionStepValue(KCPFApp.ADJOINT_AUXILIARY_VELOCITY_POTENTIAL, multiplier*current_value_aux)
