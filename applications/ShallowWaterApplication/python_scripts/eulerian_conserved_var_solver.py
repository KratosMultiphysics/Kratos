from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
# importing the Kratos Library
import KratosMultiphysics
import KratosMultiphysics.ShallowWaterApplication as Shallow

# Check that KratosMultiphysics was imported in the main script
KratosMultiphysics.CheckForPreviousImport()

## Import base class file
from shallow_water_base_solver import ShallowWaterBaseSolver

def CreateSolver(model, custom_settings):
    return EulerianConservedVarSolver(model, custom_settings)

class EulerianConservedVarSolver(ShallowWaterBaseSolver):
    def __init__(self, model, custom_settings):
        super(EulerianConservedVarSolver, self).__init__(model, custom_settings)

        # Set the element and condition names for the replace settings
        self.element_name = "EulerConsVarElement"
        self.condition_name = "Condition"
        self.min_buffer_size = 2

    def AddVariables(self):
        super(EulerianConservedVarSolver,self).AddVariables()
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.MOMENTUM)

    def AddDofs(self):
        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.MOMENTUM_X, self.main_model_part)
        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.MOMENTUM_Y, self.main_model_part)
        KratosMultiphysics.VariableUtils().AddDof(Shallow.HEIGHT, self.main_model_part)

        self.print_on_rank_zero("::[EulerianConservedVarSolver]::", "Shallow water solver DOFs added correctly.")

    def SolveSolutionStep(self):
        if self._TimeBufferIsInitialized:
            # If a node and it's neighbours are dry, set ACTIVE flag to false
            self.ShallowVariableUtils.SetDryWetState()
            # Solve equations
            is_converged = self.solver.SolveSolutionStep()
            # Compute free surface
            self.ShallowVariableUtils.ComputeFreeSurfaceElevation()
            # Compute velocity
            self.ShallowVariableUtils.ComputeVelocity()
            # If water height is negative or close to zero, reset values
            self.ShallowVariableUtils.CheckDryConservedVariables()

            return is_converged
