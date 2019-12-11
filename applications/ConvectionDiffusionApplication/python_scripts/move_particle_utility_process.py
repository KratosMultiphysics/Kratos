import KratosMultiphysics as KM
import KratosMultiphysics.ConvectionDiffusionApplication as CDA

def Factory(settings, Model):
    if not isinstance(settings, KM.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return MoveParticleUtilityProcess(Model, settings["Parameters"])


class MoveParticleUtilityProcess(KM.Process):
    def __init__(self, model, settings):

        KM.Process.__init__(self)

        default_settings = KM.Parameters("""
            {
                "model_part_name"                 : "please_specify_model_part_name",
                "use_mesh_velocity"               : false,
                "reset_boundary_conditions"       : true
            }
            """
            )

        settings.ValidateAndAssignDefaults(default_settings)

        self.model_part = model[settings["model_part_name"].GetString()]
        self.use_mesh_velocity = settings["use_mesh_velocity"].GetBool()
        self.reset_boundary_conditions = settings["reset_boundary_conditions"].GetBool()
        KM.VariableUtils().SetNonHistoricalVariableToZero(KM.MEAN_SIZE, self.model_part.Nodes)
        KM.VariableUtils().SetNonHistoricalVariableToZero(KM.MEAN_VEL_OVER_ELEM_SIZE, self.model_part.Nodes)

    def Check(self):
        settings = self.model_part.ProcessInfo[KM.CONVECTION_DIFFUSION_SETTINGS]
        if not settings.IsDefinedUnknownVariable():
            raise Exception('Unknown variable not defined in convection diffusion settings')
        if not settings.IsDefinedProjectionVariable():
            raise Exception('Projection variable not defined in convection diffusion settings')
        if not settings.IsDefinedVelocityVariable():
            raise Exception('Velocity variable not defined in convection diffusion settings')
        if not settings.IsDefinedMeshVelocityVariable():
            raise Exception('MeshVelocity variable not defined in convection diffusion settings')

    def ExecuteBeforeSolutionLoop(self):
        dimension = self.model_part.ProcessInfo[KM.DOMAIN_SIZE]
        num_of_avg_elems = 10
        num_of_avg_nodes = 10
        neighbor_search = KM.FindNodalNeighboursProcess(self.model_part, num_of_avg_elems, num_of_avg_nodes)
        neighbor_search.Execute()
        neighbor_elements_search = KM.FindElementalNeighboursProcess(self.model_part, dimension, num_of_avg_elems)
        neighbor_elements_search.Execute()

        settings = self.model_part.ProcessInfo[KM.CONVECTION_DIFFUSION_SETTINGS]
        self.unknown_var = settings.GetUnknownVariable()
        self.projection_var = settings.GetProjectionVariable()
        self.velocity_var = settings.GetVelocityVariable()
        self.mesh_velocity_var = settings.GetMeshVelocityVariable()

        max_num_of_particles = 8 * dimension
        if dimension == 2:
            self.moveparticles = CDA.MoveParticleUtilityScalarTransport2D(self.model_part, max_num_of_particles)
        else:
            self.moveparticles = CDA.MoveParticleUtilityScalarTransport3D(self.model_part, max_num_of_particles)
        self.moveparticles.MountBin()

    def ExecuteInitializeSolutionStep(self):
        if not self.use_mesh_velocity:
            KM.VariableUtils().SetVectorVar(self.mesh_velocity_var, [0.0, 0.0, 0.0], self.model_part.Nodes)
        self.moveparticles.CalculateVelOverElemSize()
        self.moveparticles.MoveParticles()
        dimension = self.model_part.ProcessInfo[KM.DOMAIN_SIZE]
        pre_minimum_num_of_particles = dimension
        self.moveparticles.PreReseed(pre_minimum_num_of_particles)
        self.moveparticles.TransferLagrangianToEulerian()
        KM.VariableUtils().CopyScalarVar(self.projection_var, self.unknown_var, self.model_part.Nodes)
        if self.reset_boundary_conditions:
            self.moveparticles.ResetBoundaryConditions()
        self.moveparticles.CopyScalarVarToPreviousTimeStep(self.unknown_var, self.model_part.Nodes)

    def ExecuteFinalizeSolutionStep(self):
        self.moveparticles.CalculateDeltaVariables()
        self.moveparticles.CorrectParticlesWithoutMovingUsingDeltaVariables()
        dimension = self.model_part.ProcessInfo[KM.DOMAIN_SIZE]
        post_minimum_num_of_particles = 2 * dimension
        self.moveparticles.PostReseed(post_minimum_num_of_particles)

