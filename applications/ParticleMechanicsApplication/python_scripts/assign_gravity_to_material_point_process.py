import KratosMultiphysics
import KratosMultiphysics.MPMApplication as KratosMPM

def Factory(settings, Model):
    if(not isinstance(settings, KratosMultiphysics.Parameters)):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return AssignGravityToMaterialPointProcess(Model, settings["Parameters"])

## All the processes python should be derived from "Process"
class AssignGravityToMaterialPointProcess(KratosMultiphysics.Process):
    def __init__(self, Model, settings ):
        KratosMultiphysics.Process.__init__(self)

        default_settings = KratosMultiphysics.Parameters("""
            {
                "mesh_id"              : 0,
                "model_part_name"      : "please_specify_model_part_name",
                "variable_name"        : "MP_VOLUME_ACCELERATION",
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

        # Detect if variable_name is MP_VOLUME_ACCELERATION
        if(settings.Has("variable_name")):
            if(settings["variable_name"].GetString() != "MP_VOLUME_ACCELERATION"):
                KratosMultiphysics.Logger.PrintInfo("Warning in apply gravity to material point", "Error in determining variable_name")
                raise Exception('The assign_gravity_to_material_point_process only accepts \"MP_VOLUME_ACCELERATION\" as variable_name.')

        settings.ValidateAndAssignDefaults(default_settings)

        # Default settings
        self.model_part = Model[settings["model_part_name"].GetString()]
        self.modulus = settings["modulus"].GetDouble()
        self.gravity_direction = settings["direction"].GetVector()
        self.gravity_acceleration = self.modulus * self.gravity_direction

    def ExecuteBeforeSolutionLoop(self):
        if not self.model_part.ProcessInfo[KratosMultiphysics.IS_RESTARTED]:
            # Assign gravity to MP after solver.Initialize() - only apply once at the beginning!
            for element in self.model_part.Elements:
                element.SetValuesOnIntegrationPoints(KratosMPM.MP_VOLUME_ACCELERATION,[self.gravity_acceleration],self.model_part.ProcessInfo)
                element.SetValuesOnIntegrationPoints(KratosMPM.MP_ACCELERATION,[self.gravity_acceleration],self.model_part.ProcessInfo)
