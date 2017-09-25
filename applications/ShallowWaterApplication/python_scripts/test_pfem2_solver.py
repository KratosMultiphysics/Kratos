from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
# importing the Kratos Library
import KratosMultiphysics
import KratosMultiphysics.ShallowWaterApplication as KratosShallow

# Check that KratosMultiphysics was imported in the main script
KratosMultiphysics.CheckForPreviousImport()

## Import base class file
import shallow_water_base_solver

def CreateSolver(model_part, custom_settings):
    return TestPfem2Solver(model_part, custom_settings)

class TestPfem2Solver(shallow_water_base_solver.ShallowWaterBaseSolver):
    def __init__(self, model_part, custom_settings):
        # Model part and solver init
        super(TestPfem2Solver,self).__init__(model_part,custom_settings)
        # Particle stage init
        super(TestPfem2Solver,self)._pfem2_init(model_part)

    def AddVariables(self):
        super(TestPfem2Solver,self).AddVariables()
        super(TestPfem2Solver,self)._AddParticleVariables()

    def AddDofs(self):
        super(TestPfem2Solver,self)._AddPrimitiveDofs()

    def Initialize(self):
        super(TestPfem2Solver,self).Initialize()

        # Creating the solution strategy for the particle stage
        self.VariableUtils = KratosMultiphysics.VariableUtils()
        #~ maximum_number_of_particles= 8*self.domain_size
        self.moveparticles = KratosShallow.MoveShallowWaterParticleUtility(self.model_part, self.settings["pfem2_settings"])
        self.moveparticles.MountBin()

    def Solve(self):
        # Move particles
        super(TestPfem2Solver,self).ExecuteParticlesUtilitiesBeforeSolve()
        # Solve equations on mesh
        #~ (self.solver).Solve()
        for node in self.model_part.Nodes:
            node.SetSolutionStepValue(KratosShallow.HEIGHT,          node.GetSolutionStepValue(KratosShallow.PROJECTED_SCALAR1) )
            node.SetSolutionStepValue(KratosMultiphysics.VELOCITY_X, node.GetSolutionStepValue(KratosShallow.PROJECTED_VECTOR1_X) )
            node.SetSolutionStepValue(KratosMultiphysics.VELOCITY_Y, node.GetSolutionStepValue(KratosShallow.PROJECTED_VECTOR1_Y) )
            node.SetSolutionStepValue(KratosMultiphysics.VELOCITY_Z, node.GetSolutionStepValue(KratosShallow.PROJECTED_VECTOR1_Z) )
        # Compute free surface
        (self.ShallowVariableUtils).ComputeFreeSurfaceElevation()
        # Update particles
        super(TestPfem2Solver,self).ExecuteParticlesUtilitiesAfetrSolve()
