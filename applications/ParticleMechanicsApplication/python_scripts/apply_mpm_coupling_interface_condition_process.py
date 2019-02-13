import KratosMultiphysics
import KratosMultiphysics.ParticleMechanicsApplication as KratosParticle
from KratosMultiphysics.ParticleMechanicsApplication.apply_mpm_particle_dirichlet_condition_process import ApplyMPMParticleDirichletConditionProcess

def Factory(settings, Model):
    if(not isinstance(settings, KratosMultiphysics.Parameters)):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return ApplyMPMCouplingInterfaceConditionProcess(Model, settings["Parameters"])

## All the processes python should be derived from "Process"
class ApplyMPMCouplingInterfaceConditionProcess(ApplyMPMParticleDirichletConditionProcess):
    def __init__(self, Model, settings ):

        default_parameters = KratosMultiphysics.Parameters( """
            {
                "model_part_name"           : "PLEASE_SPECIFY_MODEL_PART_NAME",
                "particles_per_condition"   : 0,
                "imposition_type"           : "penalty",
                "penalty_factor"            : 0,
                "constrained"               : "fixed",
                "option"                    : ""
            }  """ )

        settings.ValidateAndAssignDefaults(default_parameters)

        super(ApplyMPMCouplingInterfaceConditionProcess, self).__init__(Model, settings)

        self.model_part = Model[settings["model_part_name"].GetString()]
        KratosMultiphysics.VariableUtils().SetFlag(KratosMultiphysics.INTERFACE, True, self.model_part.Conditions)



