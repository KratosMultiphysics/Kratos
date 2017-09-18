from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
# importing the Kratos Library
import KratosMultiphysics
import KratosMultiphysics.ShallowWaterApplication as KratosShallow

# Check that KratosMultiphysics was imported in the main script
KratosMultiphysics.CheckForPreviousImport()

## Import base class file
import shallow_water_base_solver

def CreateSolver(model_part, custom_settings):
    return Pfem2PrimitiveVarSolver(model_part, custom_settings)

class Pfem2PrimitiveVarSolver(shallow_water_base_solver.ShallowWaterBaseSolver):
    def __init__(self, model_part, custom_settings):
        # Model part and solver init
        super(Pfem2PrimitiveVarSolver,self).__init__(model_part,custom_settings)
        # Particle stage init
        super(Pfem2PrimitiveVarSolver,self)._pfem2_init(model_part)

    def AddVariables(self):
        super(Pfem2PrimitiveVarSolver,self).AddVariables()
        super(Pfem2PrimitiveVarSolver,self)._AddParticleVariables()
        self.model_part.AddNodalSolutionStepVariable(KratosShallow.DELTA_VELOCITY)
        self.model_part.AddNodalSolutionStepVariable(KratosShallow.PROJECTED_VELOCITY)

    def AddDofs(self):
        super(Pfem2PrimitiveVarSolver,self)._AddPrimitiveDofs()

    def Initialize(self):
        super(Pfem2PrimitiveVarSolver,self).Initialize()

        # Creating the solution strategy for the particle stage
        self.VariableUtils = KratosMultiphysics.VariableUtils()
        maximum_number_of_particles= 8*self.domain_size
        self.moveparticles = KratosShallow.MoveShallowWaterParticleUtility(self.model_part,maximum_number_of_particles)  
        self.moveparticles.MountBin()

    def Solve(self):
        # Move particles
        super(Pfem2PrimitiveVarSolver,self).ExecuteParticlesUtilitiesBeforeSolve()
        # Solve equations on mesh
        (self.solver).Solve()
        # Compute free surface
        (self.ShallowVariableUtils).ComputeFreeSurfaceElevation()
        # Update particles
        super(Pfem2PrimitiveVarSolver,self).ExecuteParticlesUtilitiesAfetrSolve()
