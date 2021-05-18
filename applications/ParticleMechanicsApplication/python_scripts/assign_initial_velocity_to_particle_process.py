import KratosMultiphysics
import KratosMultiphysics.ParticleMechanicsApplication as KratosParticle

def Factory(settings, Model):
    if(not isinstance(settings, KratosMultiphysics.Parameters)):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return AssignInitialVelocityToParticleProcess(Model, settings["Parameters"])

## All the processes python should be derived from "Process"
class AssignInitialVelocityToParticleProcess(KratosMultiphysics.Process):
    def __init__(self, Model, settings ):
        KratosMultiphysics.Process.__init__(self)

        default_settings = KratosMultiphysics.Parameters("""
            {
                "mesh_id"              : 0,
                "model_part_name"      : "please_specify_model_part_name",
                "variable_name"        : "MP_VELOCITY",
                "modulus"              : 1.0,
                "constrained"          : true,
                "direction"            : [0.0, 0.0, 0.0],
                "local_axes"           : {}
            }
            """)

        # Trick: allow "modulus" and "direction" to be a double or a string value (otherwise the ValidateAndAssignDefaults might fail)
        if(settings.Has("modulus")):
            if(settings["modulus"].IsString()):
                default_settings["modulus"].SetString("0.0")

        if(settings.Has("direction")):
            if(settings["direction"].IsString()):
                default_settings["direction"].SetString("Automatic")

        # Detect if variable_name is MP_VELOCITY
        if(settings.Has("variable_name")):
            if(settings["variable_name"].GetString() != "MP_VELOCITY"):
                KratosMultiphysics.Logger.PrintInfo("Warning in apply velocity to particle", "Error in determining variable_name")
                raise Exception('The assign_initial_velocity_to_particle_process only accepts \"MP_VELOCITY\" as variable_name.')

        settings.ValidateAndAssignDefaults(default_settings)

        # Get updated model_part
        self.model = Model
        model_part_name = settings["model_part_name"].GetString()
        if (model_part_name.startswith('Initial_MPM_Material.')):
            model_part_name = model_part_name.replace('Initial_MPM_Material.','')
        self.mpm_material_model_part_name = "MPM_Material." + model_part_name
        # The actual initial velocity application occurs after the submodelpart is
        # transferred from the initial MPM material to the MPM material in the particle
        # generator utility. Therefore we change the prefix from initial MPM material
        # to MPM material.

        # Default settings
        self.modulus = settings["modulus"].GetDouble()
        self.velocity_direction = settings["direction"].GetVector()
        self.velocity = self.modulus * self.velocity_direction

    def ExecuteBeforeSolutionLoop(self):
        # Assign velocity to MP after solver.Initialize() - only apply once at the beginning!
        model_part = self.model[self.mpm_material_model_part_name]
        # the model part is identified here, AFTER it has been transferred to the MPM_material part!
        for element in model_part.Elements:
            element.SetValuesOnIntegrationPoints(KratosParticle.MP_VELOCITY,[self.velocity],model_part.ProcessInfo)