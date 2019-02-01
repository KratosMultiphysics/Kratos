import KratosMultiphysics
import KratosMultiphysics.ParticleMechanicsApplication as KratosParticle

def Factory(settings, Model):
    if(not isinstance(settings, KratosMultiphysics.Parameters)):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return ApplyMPMParticleDirichletConditionProcess(Model, settings["Parameters"])

## All the processes python should be derived from "Process"
class ApplyMPMParticleDirichletConditionProcess(KratosMultiphysics.Process):
    def __init__(self, Model, settings ):
        KratosMultiphysics.Process.__init__(self)

        default_parameters = KratosMultiphysics.Parameters( """
            {
                "model_part_name"           : "PLEASE_SPECIFY_MODEL_PART_NAME",
                "particles_per_condition"   : 0,
                "imposition_type"           : "penalty",
                "penalty_factor"            : 0,
                "variable_name"             : "MPC_DISPLACEMENT",
                "modulus"                   : 1.0,
                "constrained"               : true,
                "direction"                 : [0.0, 0.0, 0.0],
                "local_axes"                : {}
            }  """ )

        settings.ValidateAndAssignDefaults(default_parameters)

        self.model_part = Model[settings["model_part_name"].GetString()]
        self.particles_per_condition = settings["particles_per_condition"].GetInt()
        self.imposition_type = settings["imposition_type"].GetString()
        self.is_neumann_boundary = False

        # set type of boundary
        if (self.imposition_type == "penalty" or self.imposition_type == "Penalty"):
            self.penalty_factor = settings["penalty_factor"].GetDouble()
        else:
            err_msg =  "The requested type of Dirichlet boundary imposition: \"" + self.imposition_type + "\" is not available!\n"
            err_msg += "Available option is: \"penalty\"."
            raise Exception(err_msg)

        # get variable imposed and check
        self.variable = KratosMultiphysics.KratosGlobals.GetVariable(settings["variable_name"].GetString())
        if(type(self.variable) != KratosMultiphysics.Array1DVariable3 and type(self.variable) != KratosMultiphysics.VectorVariable):
            msg = "Error in AssignVectorVariableProcess. Variable type of variable : " + settings["variable_name"].GetString() + " is incorrect . Must be a vector or array3"
            raise Exception(msg)

        self.modulus = settings["modulus"].GetDouble()
        self.vector_direction = settings["direction"].GetVector()
        self.vector = self.modulus * self.vector_direction

        # Compute the normal on the nodes of interest -
        KratosMultiphysics.NormalCalculationUtils().CalculateOnSimplex(self.model_part, self.model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE])

        # Set Flag BOUNDARY and variables PARTICLES_PER_CONDITION
        if self.particles_per_condition >= 0:
            KratosMultiphysics.VariableUtils().SetFlag(KratosMultiphysics.BOUNDARY, True, self.model_part.Nodes)

            for condition in self.model_part.Conditions:
                condition.Set(KratosMultiphysics.BOUNDARY, True)
                condition.SetValue(KratosParticle.PARTICLES_PER_CONDITION, self.particles_per_condition)
                condition.SetValue(KratosParticle.MPC_IS_NEUMANN, self.is_neumann_boundary)
                condition.SetValue(KratosParticle.PENALTY_FACTOR, self.penalty_factor)
                condition.SetValue(self.variable, self.vector)
        else:
            err_msg = '\n::[ApplyMPMParticleDirichletConditionProcess]:: W-A-R-N-I-N-G: You have specified invalid "particles_per_condition", '
            err_msg += 'or assigned negative values. \nPlease assign: "particles_per_condition" > 0 or = 0 (for automatic value)!\n'
            raise Exception(err_msg)