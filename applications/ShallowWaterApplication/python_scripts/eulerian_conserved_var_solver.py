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
    def __init__(self, model_part, custom_settings):
        # Model part and solver init
        super(EulerianConservedVarSolver,self).__init__(model_part,custom_settings)

    def AddVariables(self):
        super(EulerianConservedVarSolver,self).AddVariables()
        self.model_part.AddNodalSolutionStepVariable(KratosMultiphysics.VELOCITY)

    def AddDofs(self):
        super(EulerianConservedVarSolver,self)._AddConservedDofs()

    def Initialize(self):
        super(EulerianConservedVarSolver,self)._InitializeMeshStage()

    def Solve(self):
        (self.solver).Solve()
