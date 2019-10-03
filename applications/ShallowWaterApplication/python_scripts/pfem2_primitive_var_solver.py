from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
# importing the Kratos Library
import KratosMultiphysics as KM
import KratosMultiphysics.ShallowWaterApplication as SW

## Import base class file
from KratosMultiphysics.ShallowWaterApplication.shallow_water_base_solver import ShallowWaterBaseSolver

def CreateSolver(model, custom_settings):
    return Pfem2PrimitiveVarSolver(model, custom_settings)

class Pfem2PrimitiveVarSolver(ShallowWaterBaseSolver):
    def __init__(self, model, settings):
        super(Pfem2PrimitiveVarSolver, self).__init__(model, settings)

        # Set the element and condition names for the replace settings
        self.element_name = "PFEM2ReducedSWE"
        self.condition_name = "Condition"
        self.min_buffer_size = 2

        # Pfem2 settings
        domain_size = self.main_model_part.ProcessInfo[KM.DOMAIN_SIZE]
        self.settings.AddValue("pfem2_settings", KM.Parameters("""{}"""))
        self.settings["pfem2_settings"].AddEmptyValue("convection_scalar_variable").SetString("HEIGHT")
        self.settings["pfem2_settings"].AddEmptyValue("convection_vector_variable").SetString("VELOCITY")
        self.settings["pfem2_settings"].AddEmptyValue("maximum_number_of_particles").SetInt(8*domain_size)

        self.print_particles = False # By default we don't print the particles, it is expensive
        if self.print_particles:
            self.lagrangian_model_part = model.CreateModelPart("pfem2_particles")
            self.filter_factor = 1

    def AddVariables(self):
        super(Pfem2PrimitiveVarSolver, self).AddVariables()
        # Variables to project unknown and update particles
        self.main_model_part.AddNodalSolutionStepVariable(SW.DELTA_SCALAR1)
        self.main_model_part.AddNodalSolutionStepVariable(SW.DELTA_VECTOR1)
        # Specific variables to convect particles
        self.main_model_part.AddNodalSolutionStepVariable(KM.YP)
        self.main_model_part.AddNodalSolutionStepVariable(SW.MEAN_SIZE)

    def AddDofs(self):
        KM.VariableUtils().AddDof(KM.VELOCITY_X, self.main_model_part)
        KM.VariableUtils().AddDof(KM.VELOCITY_Y, self.main_model_part)
        KM.VariableUtils().AddDof(SW.HEIGHT, self.main_model_part)

        KM.Logger.PrintInfo("::[Pfem2PrimitiveVarSolver]::", "Shallow water solver DOFs added correctly.")

    def Initialize(self):
        super(Pfem2PrimitiveVarSolver, self).Initialize()

        # Initializing the neighbour search
        domain_size = self.main_model_part.ProcessInfo[KM.DOMAIN_SIZE]
        number_of_avg_elems = 10
        number_of_avg_nodes = 10
        self.neighbour_search = KM.FindNodalNeighboursProcess(self.main_model_part, number_of_avg_elems, number_of_avg_nodes)
        self.neighbour_search.Execute()
        self.neighbour_elements_search = KM.FindElementalNeighboursProcess(self.main_model_part, domain_size, number_of_avg_elems)
        self.neighbour_elements_search.Execute()

        # Creating the solution strategy for the particle stage
        self.moveparticles = SW.MoveShallowWaterParticleUtility(self.main_model_part, self.settings["pfem2_settings"])
        self.moveparticles.MountBin()
        if self.print_particles:
            self.moveparticles.ExecuteParticlesPrintingTool(self.lagrangian_model_part, self.filter_factor)
        KM.Logger.PrintInfo("::[Pfem2PrimitiveVarSolver]::", "Pfem2 stage initialization finished")

    def InitializeSolutionStep(self):
        if self._TimeBufferIsInitialized():
            # Move particles
            self.moveparticles.CalculateVelOverElemSize()
            self.moveparticles.MoveParticles()
            # Reseed empty elements
            pre_minimum_number_of_particles = self.main_model_part.ProcessInfo[KM.DOMAIN_SIZE]
            self.moveparticles.PreReseed(pre_minimum_number_of_particles)
            # Project info to mesh
            self.moveparticles.TransferLagrangianToEulerian()
            self.moveparticles.ResetBoundaryConditions()
            # Initialize mesh solution step
            self.solver.InitializeSolutionStep()

    def FinalizeSolutionStep(self):
        if self._TimeBufferIsInitialized():
            # Finalize mesh solution step
            self.solver.FinalizeSolutionStep()

            # Update particles
            self.moveparticles.CalculateDeltaVariables()
            self.moveparticles.CorrectParticlesWithoutMovingUsingDeltaVariables()
            # Reseed empty elements
            post_minimum_number_of_particles = self.main_model_part.ProcessInfo[KM.DOMAIN_SIZE]*2
            self.moveparticles.PostReseed(post_minimum_number_of_particles)

            # Print particles if needed
            if self.print_particles:
                self.lagrangian_model_part.ProcessInfo[KratosMultiphysics.STEP] = self.main_model_part.ProcessInfo[KratosMultiphysics.STEP]
                self.lagrangian_model_part.ProcessInfo[KratosMultiphysics.TIME] = self.main_model_part.ProcessInfo[KratosMultiphysics.TIME]
                self.moveparticles.ExecuteParticlesPrintingTool(self.lagrangian_model_part, self.filter_factor)