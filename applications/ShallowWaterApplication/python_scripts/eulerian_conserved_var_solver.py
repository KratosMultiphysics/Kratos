from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
# importing the Kratos Library
import KratosMultiphysics
import KratosMultiphysics.ShallowWaterApplication as KratosShallow

# Check that KratosMultiphysics was imported in the main script
KratosMultiphysics.CheckForPreviousImport()

## Import base class file
import shallow_water_base_solver

def CreateSolver(model_part, custom_settings):
    return EulerianConservedVarSolver(model_part, custom_settings)

class EulerianConservedVarSolver(shallow_water_base_solver.ShallowWaterBaseSolver):
    def AddVariables(self):
        super(EulerianConservedVarSolver,self).AddVariables()
        self.model_part.AddNodalSolutionStepVariable(KratosMultiphysics.MOMENTUM)

    def AddDofs(self):
        super(EulerianConservedVarSolver,self)._AddConservedDofs()

    def Solve(self):
        # If a node and it's neighbours are dry, set ACTIVE flag to false
        (self.ShallowVariableUtils).SetDryWetState()
        # Solve equations
        (self.solver).Solve()
        # Compute free surface
        (self.ShallowVariableUtils).ComputeFreeSurfaceElevation()
        # Compute velocity
        (self.ShallowVariableUtils).ComputeVelocity()
        # If water height is negative or close to zero, reset values
        (self.ShallowVariableUtils).CheckDryConservedVariables()
