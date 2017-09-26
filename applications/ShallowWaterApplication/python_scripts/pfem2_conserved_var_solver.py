from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
# importing the Kratos Library
import KratosMultiphysics
import KratosMultiphysics.ShallowWaterApplication as KratosShallow

# Check that KratosMultiphysics was imported in the main script
KratosMultiphysics.CheckForPreviousImport()

## Import base class file
import shallow_water_base_solver
import pure_convection_solver

def CreateSolver(model_part, custom_settings):
    return Pfem2ConservedVarSolver(model_part, custom_settings)

class Pfem2ConservedVarSolver(pure_convection_solver.PureConvectionSolver,shallow_water_base_solver.ShallowWaterBaseSolver):
    def __init__(self, model_part, custom_settings):
        super(Pfem2ConservedVarSolver,self).__init__(model_part,custom_settings)

    def AddVariables(self):
        super(Pfem2ConservedVarSolver,self).AddVariables()
        self.model_part.AddNodalSolutionStepVariable(KratosMultiphysics.MOMENTUM)

    def AddDofs(self):
        super(Pfem2ConservedVarSolver,self)._AddConservedDofs()

    def Initialize(self):
        shallow_water_base_solver.ShallowWaterBaseSolver.Initialize(self)
        super(Pfem2ConservedVarSolver,self).Initialize()

    def Solve(self):
        # Move particles
        super(Pfem2ConservedVarSolver,self)._ExecuteParticlesUtilityBeforeSolve()
        # Solve equations on mesh
        (self.solver).Solve()
        # Compute free surface
        (self.ShallowVariableUtils).ComputeFreeSurfaceElevation()
        # Compute velocity
        (self.ShallowVariableUtils).ComputeVelocity()
        # Update particles
        super(Pfem2ConservedVarSolver,self)._ExecuteParticlesUtilityAfterSolve()
