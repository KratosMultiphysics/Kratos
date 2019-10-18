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

        # self.response_function_settings = custom_settings["response_function_settings"].Clone()
        # self.sensitivity_settings = custom_settings["sensitivity_settings"].Clone()
        # custom_settings.RemoveValue("response_function_settings")
        # custom_settings.RemoveValue("sensitivity_settings")
        # Construct the base solver.
        super(PotentialFlowStochasticAdjointSolver, self).__init__(model, custom_settings)

        KratosMultiphysics.Logger.PrintInfo("::[PotentialFlowStochasticAdjointSolver]:: ", "Construction finished")


    def FinalizeSolutionStep(self):

        cvar_beta = self.response_function_settings["cvar_beta"].GetDouble()
        cvar_t = self.response_function_settings["cvar_t"].GetDouble()

        ################### SWITCHING TO NEGATIVE ##########################################
        objective_function = - self.response_function.CalculateValue(self.main_model_part)
        ################### SWITCHING TO NEGATIVE ##########################################

        if not (objective_function < cvar_t):
            multiplier = 1.0 / (1.0-cvar_beta)
        else:
            multiplier = 0.0
        print("CVAR_T VALUE:", cvar_t)

        print("LIFT:", objective_function)

        print("MULPLIER:", multiplier)
        for node in self.main_model_part.Nodes:
            current_value = node.GetSolutionStepValue(KCPFApp.ADJOINT_VELOCITY_POTENTIAL)
            current_value_aux = node.GetSolutionStepValue(KCPFApp.ADJOINT_AUXILIARY_VELOCITY_POTENTIAL)
            node.SetSolutionStepValue(KCPFApp.ADJOINT_VELOCITY_POTENTIAL, multiplier*current_value)
            node.SetSolutionStepValue(KCPFApp.ADJOINT_AUXILIARY_VELOCITY_POTENTIAL, multiplier*current_value_aux)

        super(PotentialFlowStochasticAdjointSolver, self).FinalizeSolutionStep()
