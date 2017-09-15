from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
# importing the Kratos Library
import KratosMultiphysics
import KratosMultiphysics.ShallowWaterApplication as KratosShallow

# Check that KratosMultiphysics was imported in the main script
KratosMultiphysics.CheckForPreviousImport()

## Import base class file
import shallow_water_base_solver

def CreateSolver(model_part, custom_settings):
    return LagrangianConservedVarSolver(model_part, custom_settings)

class LagrangianConservedVarSolver(shallow_water_base_solver.ShallowWaterBaseSolver):
    def __init__(self, model_part, custom_settings):
        # Model part and solver init
        super(LagrangianConservedVarSolver,self).__init__(model_part,custom_settings)
        # Particle stage init
        super(LagrangianConservedVarSolver,self)._pfem2_init(model_part)

    def AddVariables(self):
        super(LagrangianConservedVarSolver,self).AddVariables()
        super(LagrangianConservedVarSolver,self)._AddParticleVariables()
        self.model_part.AddNodalSolutionStepVariable(KratosMultiphysics.MOMENTUM)
        self.model_part.AddNodalSolutionStepVariable(KratosShallow.DELTA_MOMENTUM)
        self.model_part.AddNodalSolutionStepVariable(KratosShallow.PROJECTED_MOMENTUM)

    def AddDofs(self):
        super(LagrangianConservedVarSolver,self)._AddConservedDofs()

    def Initialize(self):
        super(LagrangianConservedVarSolver,self).Initialize()

        # Creating the solution strategy for the particle stage
        self.VariableUtils = KratosMultiphysics.VariableUtils()
        maximum_number_of_particles= 8*self.domain_size
        self.moveparticles = KratosShallow.MoveShallowWaterParticleHUUtility(self.model_part,maximum_number_of_particles)  
        self.moveparticles.MountBin()

    def Solve(self):
        # Move particles
        super(LagrangianConservedVarSolver,self).ExecuteParticlesUtilitiesBeforeSolve()
        # Solve equations on mesh
        (self.solver).Solve()
        # Update particles
        super(LagrangianConservedVarSolver,self).ExecuteParticlesUtilitiesAfetrSolve()

