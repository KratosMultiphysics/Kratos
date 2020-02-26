import KratosMultiphysics as KM
import KratosMultiphysics.FluidTransportApplication as FTA

from KratosMultiphysics.ConvectionDiffusionApplication.move_particle_utility_process import MoveParticleUtilityProcess

def Factory(settings, Model):
    if not isinstance(settings, KM.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return PFEM2FluidTransportProcess(Model, settings["Parameters"])


class PFEM2FluidTransportProcess(MoveParticleUtilityProcess):
    def __init__(self, model, settings):

        default_settings = KM.Parameters("""
            {
                "model_part_name"                 : "please_specify_model_part_name",
                "use_mesh_velocity"               : false,
                "reset_boundary_conditions"       : true,
                "crank_nicolson_theta"            : 0.5
            }
            """
            )

        settings.ValidateAndAssignDefaults(default_settings)

        self.model_part = model[settings["model_part_name"].GetString()]
        self.use_mesh_velocity = settings["use_mesh_velocity"].GetBool()
        self.reset_boundary_conditions = settings["reset_boundary_conditions"].GetBool()
        self.crank_nicolson_theta = settings["crank_nicolson_theta"].GetDouble()

    def ExecuteInitializeSolutionStep(self):
    #     # Move the particles until the half step predictor
    #     time_step = self.model_part.ProcessInfo.GetValue(KM.DELTA_TIME)
    #     self.model_part.SetValue(KM.DELTA_TIME, time_step * self.crank_nicolson_theta)
        super(PFEM2FluidTransportProcess, self).ExecuteInitializeSolutionStep()
    #     self.model_part.SetValue(KM.DELTA_TIME, time_step)

        # Copy Unknown var to TEMPERATURE on the previous time step (convected) before solving the system
        KM.VariableUtils().CopyModelPartNodalVar(
            FTA.PHI_THETA,
            KM.TEMPERATURE,
            self.model_part,
            self.model_part,
            1)

    def ExecuteFinalizeSolutionStep(self):
        # Finish the movement of the particles to the final time step
        # time_step = self.model_part.ProcessInfo.GetValue(KM.DELTA_TIME)
        # self.model_part.ProcessInfo.SetValue(KM.DELTA_TIME, time_step * (1 - self.crank_nicolson_theta))
        # self.moveparticles.MoveParticles()
        # self.model_part.ProcessInfo.SetValue(KM.DELTA_TIME, time_step)

        # Copy again the TEMPERATURE to Unknown var before updating the particles
        KM.VariableUtils().CopyModelPartNodalVar(
            KM.TEMPERATURE,
            FTA.PHI_THETA,
            self.model_part,
            self.model_part,
            0)

        # update the particles with the diffusion from the mesh stage
        super(PFEM2FluidTransportProcess, self).ExecuteFinalizeSolutionStep()
