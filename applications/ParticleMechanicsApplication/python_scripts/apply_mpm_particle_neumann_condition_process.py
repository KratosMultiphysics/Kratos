import KratosMultiphysics
import KratosMultiphysics.ParticleMechanicsApplication as KratosParticle

def Factory(settings, Model):
    if(not isinstance(settings, KratosMultiphysics.Parameters)):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return ApplyMPMParticleNeumannConditionProcess(Model, settings["Parameters"])

## All the processes python should be derived from "Process"
class ApplyMPMParticleNeumannConditionProcess(KratosMultiphysics.Process):
    def __init__(self, Model, settings ):
        KratosMultiphysics.Process.__init__(self)

        default_parameters = KratosMultiphysics.Parameters( """
            {
                "model_part_name"           : "PLEASE_SPECIFY_MODEL_PART_NAME",
                "particles_per_condition"   : 0,
                "variable_name"             : "PLEASE_SPECIFY_LOADING_CONDITION",
                "modulus"                   : 1.0,
                "constrained"               : "fixed",
                "direction"                 : [0.0, 0.0, 0.0],
                "option"                    : "",
                "local_axes"                : {}
            }  """ )

        settings.ValidateAndAssignDefaults(default_parameters)

        self.model_part = Model[settings["model_part_name"].GetString()]
        self.particles_per_condition = settings["particles_per_condition"].GetInt()
        self.is_neumann_boundary = True
        self.option = settings["option"].GetString()

        # check constraint
        self.constrained = settings["constrained"].GetString()
        if (self.constrained == "fixed"):
            self.normal_following_load = False
        elif (self.constrained == "normal"):
            self.normal_following_load = True
        else:
            err_msg =  "The requested type of constrain: \"" + self.constrained + "\" is not available!\n"
            err_msg += "Available options are: \"fixed\" and \"normal\"."
            raise Exception(err_msg)

        # get variable imposed and check
        variable_name = settings["variable_name"].GetString()
        variable_name_list = ["POINT_LOAD","LINE_LOAD","SURFACE_LOAD"]
        if(variable_name in variable_name_list):
            self.variable = KratosMultiphysics.KratosGlobals.GetVariable(variable_name)
        else:
            err_msg =  "The given variable \"" + variable_name + "\" is not available to be imposed with this process.\n"
            err_msg += "Available options are: " + ", ".join(variable_name_list)
            raise Exception(err_msg)

        self.modulus = settings["modulus"].GetDouble()
        self.vector_direction = settings["direction"].GetVector()
        self.vector = self.modulus * self.vector_direction

        # Compute the normal on the nodes of interest -
        KratosMultiphysics.NormalCalculationUtils().CalculateOnSimplex(self.model_part, self.model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE])
        self.modified_normal = False
        if self.option == "flip_normal":
            self.modified_normal = True

        # Set Flag BOUNDARY and variables PARTICLES_PER_CONDITION
        if self.particles_per_condition >= 0:
            KratosMultiphysics.VariableUtils().SetFlag(KratosMultiphysics.BOUNDARY, True, self.model_part.Nodes)

            for condition in self.model_part.Conditions:
                condition.Set(KratosMultiphysics.BOUNDARY, True)
                condition.Set(KratosMultiphysics.MARKER, self.normal_following_load)
                condition.Set(KratosMultiphysics.MODIFIED, self.modified_normal)
                condition.SetValue(KratosParticle.PARTICLES_PER_CONDITION, self.particles_per_condition)
                condition.SetValue(KratosParticle.MPC_IS_NEUMANN, self.is_neumann_boundary)
                condition.SetValue(self.variable, self.vector)
        else:
            err_msg = '\n::[ApplyMPMParticleNeumannConditionProcess]:: W-A-R-N-I-N-G: You have specified invalid "particles_per_condition", '
            err_msg += 'or assigned negative values. \nPlease assign: "particles_per_condition" > 0 or = 0 (for automatic value)!\n'
            raise Exception(err_msg)