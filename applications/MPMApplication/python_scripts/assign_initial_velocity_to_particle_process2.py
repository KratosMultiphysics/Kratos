import KratosMultiphysics
import KratosMultiphysics.MPMApplication as KratosMPM
import math

def Factory(settings, Model):
    if(not isinstance(settings, KratosMultiphysics.Parameters)):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return AssignInitialVelocityToParticleProcess2(Model, settings["Parameters"])

## All the processes python should be derived from "Process"
class AssignInitialVelocityToParticleProcess2(KratosMultiphysics.Process):
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

        # Default settings
        self.model_part = Model[settings["model_part_name"].GetString()]
        self.modulus = settings["modulus"].GetDouble()
        self.velocity_direction = settings["direction"].GetVector()
        self.velocity = self.modulus * self.velocity_direction
        gravity = 9.81
        self.height = 1.0
        frequency = 0.89
        self.wave_number = math.pi
        self.constant = gravity * self.modulus / ( 2.0 * frequency * math.cosh( self.wave_number * self.height ) ) * self.wave_number

    def ExecuteBeforeSolutionLoop(self):
        # Assign velocity to MP after solver.Initialize() - only apply once at the beginning!
        for element in self.model_part.Elements:
            mp_coord = element.CalculateOnIntegrationPoints(KratosMPM.MP_COORD,self.model_part.ProcessInfo)[0]
            self.velocity[0] = -1.0 * self.constant * math.sin( self.wave_number * mp_coord[0] ) * math.cosh( self.wave_number * ( mp_coord[1] ) )
            self.velocity[1] = self.constant * math.cos( self.wave_number * mp_coord[0] ) * math.sinh( self.wave_number * ( mp_coord[1] ) )
            self.velocity[2] = 0.0
            element.SetValuesOnIntegrationPoints(KratosMPM.MP_VELOCITY,[self.velocity],self.model_part.ProcessInfo)
