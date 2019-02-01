import KratosMultiphysics
import KratosMultiphysics.ParticleMechanicsApplication as KratosParticle

def Factory(settings, Model):
    if(not isinstance(settings, KratosMultiphysics.Parameters)):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return ApplyMPMParticleConditionProcess(Model, settings["Parameters"])

## All the processes python should be derived from "Process"
class ApplyMPMParticleConditionProcess(KratosMultiphysics.Process):
    def __init__(self, Model, settings ):
        KratosMultiphysics.Process.__init__(self)

        default_parameters = KratosMultiphysics.Parameters( """
            {
                "model_part_name":"PLEASE_CHOOSE_MODEL_PART_NAME",
                "particles_per_condition" : 0,
                "boundary_type" : ""
            }  """ )

        settings.ValidateAndAssignDefaults(default_parameters)

        self.model_part = Model[settings["model_part_name"].GetString()]
        self.particles_per_condition = settings["particles_per_condition"].GetInt()
        self.boundary_type = settings["boundary_type"].GetString()

        # set type of boundary
        if (self.boundary_type == "neumann" or self.boundary_type == "Neumann"):
            self.is_neumann_boundary = True
        elif (self.boundary_type == "dirichlet" or self.boundary_type == "Dirichlet"):
            self.is_neumann_boundary = False
        else:
            err_msg =  "The requested type of boundary \"" + self.boundary_type + "\" is not available!\n"
            err_msg += "Available options are: \"dirichlet\" and \"neumann\"."
            raise Exception(err_msg)

        # Compute the normal on the nodes of interest -
        # Note that the model part employed here is supposed to only have slip "conditions"
        KratosMultiphysics.NormalCalculationUtils().CalculateOnSimplex(self.model_part, self.model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE])

        # Set Flag BOUNDARY and variables PARTICLES_PER_CONDITION
        if self.particles_per_condition >= 0:
            KratosMultiphysics.VariableUtils().SetFlag(KratosMultiphysics.BOUNDARY, True, self.model_part.Nodes)

            for condition in self.model_part.Conditions:
                condition.Set(KratosMultiphysics.BOUNDARY, True)
                condition.SetValue(KratosParticle.PARTICLES_PER_CONDITION, self.particles_per_condition)
                condition.SetValue(KratosParticle.MPC_IS_NEUMANN, self.is_neumann_boundary)
        else:
            err_msg = '\n::[ApplyMPMParticleConditionProcess]:: W-A-R-N-I-N-G: You have specified invalid "particles_per_condition", '
            err_msg += 'or assigned negative values. \nPlease assign: "particles_per_condition" > 0 or = 0 (for automatic value)!\n'
            raise Exception(err_msg)