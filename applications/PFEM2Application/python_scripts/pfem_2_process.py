import KratosMultiphysics as KM
import KratosMultiphysics.PFEM2Application as PFEM2

def Factory(settings, Model):
    if not isinstance(settings, KM.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return PFEM2Process(Model, settings["Parameters"])


class PFEM2Process(KM.Process):
    def __init__(self, model, settings):

        KM.Process.__init__(self)

        default_settings = KM.Parameters("""
            {
                "model_part_name"                 : "please_specify_model_part_name",
                "echo_level"                      : 0,
                "use_mesh_velocity"               : false,
                "discriminate_streamlines"        : true,
                "reset_boundary_conditions"       : false,
                "fully_reset_boundary_conditions" : true,
                "print_particles"                 : false,
                "printing_model_part_name"        : "pfem2_particles",
                "printing_filter_factor"          : 5
            }
            """
            )

        settings.ValidateAndAssignDefaults(default_settings)

        self.model_part = model[settings["model_part_name"].GetString()]
        self.echo_level = settings["echo_level"].GetDouble()
        self.use_mesh_velocity = settings["use_mesh_velocity"].GetBool()
        self.discriminate_streamlines = settings["discriminate_streamlines"].GetBool()
        self.reset_boundary_conditions = settings["reset_boundary_conditions"].GetBool()
        self.full_reset_boundary_conditions = settings["fully_reset_boundary_conditions"].GetBool()

        self.print_particles = settings["print_particles"].GetBool()
        if self.print_particles:
            lagrangian_model_part_name = settings["printing_model_part_name"].GetString()
            self.lagrangian_model_part = model.CreateModelPart(lagrangian_model_part_name)
            self.filter_factor = settings["printing_filter_factor"].GetInt()

    def ExecuteInitialize(self):
        pass

    def ExecuteBeforeSolutionLoop(self):
        dimension = self.model_part.ProcessInfo[KM.DOMAIN_SIZE]
        num_of_avg_elems = 10
        num_of_avg_nodes = 10
        neighbor_search = KM.FindNodalNeighboursProcess(self.model_part)
        neighbor_search.Execute()
        neighbor_elements_search = KM.FindElementalNeighboursProcess(self.model_part, dimension, num_of_avg_elems)
        neighbor_elements_search.Execute()

        max_num_of_particles = 8 * dimension
        if dimension == 2:
            self.moveparticles = PFEM2.MoveParticleUtilityPFEM22D(self.model_part, max_num_of_particles)
        else:
            self.moveparticles = PFEM2.MoveParticleUtilityPFEM23D(self.model_part, max_num_of_particles)
        self.moveparticles.MountBin()

        if self.print_particles:
            self.moveparticles.ExecuteParticlesPritingTool(self.lagrangian_model_part, self.filter_factor)

        self.initial_water_volume = self._get_water_volume_utility().Calculate()

    def ExecuteInitializeSolutionStep(self):
        if self.use_mesh_velocity == False:
            KM.VariableUtils().SetVectorVar(KM.MESH_VELOCITY, [0.0, 0.0, 0.0], self.model_part.Nodes)
        self.moveparticles.CalculateVelOverElemSize()
        self.moveparticles.MoveParticles(self.discriminate_streamlines)
        dimension = self.model_part.ProcessInfo[KM.DOMAIN_SIZE]
        pre_minimum_num_of_particles = dimension
        self.moveparticles.PreReseed(pre_minimum_num_of_particles)
        self.moveparticles.TransferLagrangianToEulerian()
        KM.VariableUtils().CopyVectorVar(PFEM2.PROJECTED_VELOCITY, KM.VELOCITY, self.model_part.Nodes)
        if self.reset_boundary_conditions:
            self.moveparticles.ResetBoundaryConditions(self.full_reset_boundary_conditions)
        self.moveparticles.CopyVectorVarToPreviousTimeStep(KM.VELOCITY, self.model_part.Nodes)

    def ExecuteFinalizeSolutionStep(self):
        self.moveparticles.CalculateDeltaVelocity()
        self.moveparticles.AccelerateParticlesWithoutMovingUsingDeltaVelocity()
        dimension = self.model_part.ProcessInfo[KM.DOMAIN_SIZE]
        post_minimum_num_of_particles = 2 * dimension
        mass_correction_factor = self.compute_mass_correction_factor()
        self.moveparticles.PostReseed(post_minimum_num_of_particles, mass_correction_factor)

        if self.print_particles:
            self.lagrangian_model_part.ProcessInfo[KM.STEP] = self.model_part.ProcessInfo[KM.STEP]
            self.lagrangian_model_part.ProcessInfo[KM.TIME] = self.model_part.ProcessInfo[KM.TIME]
            self.moveparticles.ExecuteParticlesPritingTool(self.lagrangian_model_part, self.filter_factor)

    def compute_mass_correction_factor(self):
        water_volume = self._get_water_volume_utility().Calculate()
        water_fraction = water_volume / self.initial_water_volume if self.initial_water_volume else 1.0
        factor = 1.0
        mass_correction_factor = (1.0 - water_fraction) * factor
        if self.echo_level > 0:
            KM.Logger.PrintInfo("PFEM2Process", "Mass correction factor : ", mass_correction_factor)
        return mass_correction_factor

    def _get_water_volume_utility(self):
        if not hasattr(self, '_water_volume_utility'):
            self._water_volume_utility = self._create_water_volume_utility()
        return self._water_volume_utility

    def _create_water_volume_utility(self):
        water_volume_utility = None
        dimension = self.model_part.ProcessInfo[KM.DOMAIN_SIZE]
        if dimension == 2:
            water_volume_utility = PFEM2.CalculateWaterFraction2D(self.model_part)
        else:
            water_volume_utility = PFEM2.CalculateWaterFraction3D(self.model_part)
        return water_volume_utility
