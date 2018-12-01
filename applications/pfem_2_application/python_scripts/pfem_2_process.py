import KratosMultiphysics as Kratos
import KratosMultiphysics.KratosPFEM2Application as Pfem2

def Factory(settings, Model):
    if(type(settings) != KratosMultiphysics.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return PFEM2Process(Model, settings["Parameters"])


class PFEM2Process(KratosMultiphysics.Process):
    def __init__(self, Model, settings):

        KratosMultiphysics.Process.__init__(self)

        default_settings = KratosMultiphysics.Parameters("""
            {
                "model_part_name"      : "please_specify_model_part_name",
                "variables_to_convect" : [""],
                "projection_variables" : [""],
                "maximum_number_of_particles" : 16
            }
            """
            )

        settings.ValidateAndAssignDefaults(default_settings)

        self.use_mesh_velocity = False
        self.model_part = Model[settings["model_part_name"].GetString()]
        self.max_num_of_particles = settings["maximum_number_of_particles"].GetInt()
        self.discriminate_streamlines = True
        self.pre_minimum_num_of_particles = self.model_part.ProcessInfo[Kratos.DOMAIN_SIZE]
        self.post_minimum_num_of_particles = 2 * self.model_part.ProcessInfo[Kratos.DOMAIN_SIZE]
        self.full_reset_boundary_conditions = True
        self.mass_correction_factor = 0.0

    def ExecuteInitialize(self):
        num_of_avg_elems = 10
        num_of_avg_nodes = 10
        neighbour_search = Kratos.FindNodalNeighboursProcess(self.model_part, num_of_avg_elems, num_of_avg_nodes)
        neighbour_search.Execute()
        neighbour_elements_search = Kratos.FindElementalNeighboursProcess(self.model_part, self.dimension, num_of_avg_elems)
        neighbour_elements_search.Execute()

        if self.model_part.ProcessInfo[Kratos.DOMAIN_SIZE] == 2:
            self.moveparticles = Pfem2.MoveParticleUtilityPFEM22D(self.main_model_part, self.max_num_of_particles)
        else:
            self.moveparticles = Pfem2.MoveParticleUtilityPFEM22D(self.main_model_part, self.max_num_of_particles)
        self.moveparticles.MountBin()

    def ExecuteBeforeSolutionLoop(self):
        pass

    def ExecuteInitializeSolutionStep(self):
        if self.use_mesh_velocity == False:
            Kratos.VariableUtils().SetVectorVar(Kratos.MESH_VELOCITY, 0.0, self.model_part.Nodes)
        self.moveparticles.CalculateVelOverElemSize()
        self.moveparticles.MoveParticles(self.discriminate_streamlines)
        self.moveparticles.PreReseed(self.pre_minimum_num_of_particles)
        self.moveparticles.TransferLagrangianToEulerian()
        self.moveparticles.ResetBoundaryConditions(self.full_reset_boundary_conditions)
        Kratos.VariableUtils().CopyVectorVar(Pfem2.PROJECTED_VELOCITY, Kratos.VELOCITY, self.model_part.Nodes)
        self.moveparticles.CopyVectorVarToPreviousTimeStep(Kratos.VELOCITY, self.model_part)

    def ExecuteFinalizeSolutionStep(self):
        self.moveparticles.CalculateDeltaVelocity()
        self.moveparticles.AccelerateParticlesWithoutMovingUsingDeltaVelocity()
        self.moveparticles.PostReseed(self.post_minimum_num_of_particles, self.mass_correction_factor)
