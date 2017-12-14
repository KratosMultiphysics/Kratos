from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
# importing the Kratos Library
import KratosMultiphysics
import KratosMultiphysics.ShallowWaterApplication as KratosShallow

# Check that KratosMultiphysics was imported in the main script
KratosMultiphysics.CheckForPreviousImport()

## Import base class file
import shallow_water_base_solver

def CreateSolver(model_part, custom_settings):
    return PureConvectionSolver(model_part, custom_settings)

class PureConvectionSolver(shallow_water_base_solver.ShallowWaterBaseSolver):
    def __init__(self, model_part, custom_settings):
        # Model part and solver init
        super(PureConvectionSolver,self).__init__(model_part,custom_settings)
        # For the pfem2
        number_of_avg_elems = 10
        number_of_avg_nodes = 10
        self.neighbour_search = KratosMultiphysics.FindNodalNeighboursProcess(model_part,number_of_avg_elems,number_of_avg_nodes)
        (self.neighbour_search).Execute()
        self.neighbour_elements_search= KratosMultiphysics.FindElementalNeighboursProcess(model_part,self.domain_size,number_of_avg_elems)
        (self.neighbour_elements_search).Execute()

    def AddVariables(self):
        super(PureConvectionSolver,self).AddVariables()
        # Variables to project unknown and update particles
        self.model_part.AddNodalSolutionStepVariable(KratosShallow.DELTA_SCALAR1)
        self.model_part.AddNodalSolutionStepVariable(KratosShallow.PROJECTED_SCALAR1)
        self.model_part.AddNodalSolutionStepVariable(KratosShallow.DELTA_VECTOR1)
        self.model_part.AddNodalSolutionStepVariable(KratosShallow.PROJECTED_VECTOR1)
        # Specific variables to convect particles
        self.model_part.AddNodalSolutionStepVariable(KratosMultiphysics.YP);
        self.model_part.AddNodalSolutionStepVariable(KratosShallow.MEAN_SIZE);

    def AddDofs(self):
        super(PureConvectionSolver,self)._AddPrimitiveDofs()

    def Initialize(self):
        # Creating the solution strategy for the particle stage
        self.VariableUtils = KratosMultiphysics.VariableUtils()
        #~ maximum_number_of_particles= 8*self.domain_size
        self.moveparticles = KratosShallow.MoveShallowWaterParticleUtility(self.model_part, self.settings["pfem2_settings"])
        self.moveparticles.MountBin()
        print("Pfem2 utility initialization finished")

    def Solve(self):
        # Move particles
        self._ExecuteParticlesUtilityBeforeSolve()
        # Solve equations on mesh: no solver, so, copy variables
        (self.VariableUtils).CopyScalarVar(KratosShallow.PROJECTED_SCALAR1,KratosShallow.HEIGHT,self.model_part.Nodes)
        (self.VariableUtils).CopyVectorVar(KratosShallow.PROJECTED_VECTOR1,KratosMultiphysics.VELOCITY,self.model_part.Nodes)
        # Compute free surface
        (self.ShallowVariableUtils).ComputeFreeSurfaceElevation()
        # Update particles
        self._ExecuteParticlesUtilityAfterSolve()

    def _ExecuteParticlesUtilityBeforeSolve(self):
        # Move particles
        (self.moveparticles).CalculateVelOverElemSize();
        (self.moveparticles).MoveParticles();
        pre_minimum_number_of_particles=self.domain_size;
        (self.moveparticles).PreReseed(pre_minimum_number_of_particles);
        (self.moveparticles).TransferLagrangianToEulerian();
        #(self.VariableUtils).CopyScalarVar(KratosShallow.PROJECTED_SCALAR1,KratosShallow.HEIGHT,self.model_part.Nodes)
        #(self.VariableUtils).CopyVectorVar(KratosShallow.PROJECTED_VECTOR1,KratosMultiphysics.VELOCITY,self.model_part.Nodes)
        (self.moveparticles).ResetBoundaryConditions()
        #(self.moveparticles).CopyScalarVarToPreviousTimeStep(KratosShallow.PROJECTED_SCALAR1,self.model_part.Nodes)
        #(self.moveparticles).CopyVectorVarToPreviousTimeStep(KratosShallow.PROJECTED_VECTOR1,self.model_part.Nodes)

    def _ExecuteParticlesUtilityAfterSolve(self):
        # Update particles
        (self.moveparticles).CalculateDeltaVariables();
        (self.moveparticles).CorrectParticlesWithoutMovingUsingDeltaVariables();
        post_minimum_number_of_particles=self.domain_size*2;
        (self.moveparticles).PostReseed(post_minimum_number_of_particles);
