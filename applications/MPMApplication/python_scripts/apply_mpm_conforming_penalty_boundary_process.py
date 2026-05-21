import KratosMultiphysics
import KratosMultiphysics.MPMApplication as KratosMPM

def Factory(settings, Model):
    if(not isinstance(settings, KratosMultiphysics.Parameters)):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return ApplyMPMConformingPenaltyBoundaryProcess(Model, settings["Parameters"])

class ApplyMPMConformingPenaltyBoundaryProcess(KratosMultiphysics.Process):
    def __init__(self, Model, settings ):
        KratosMultiphysics.Process.__init__(self)

        default_parameters = KratosMultiphysics.Parameters( """
            {
                "model_part_name"                : "PLEASE_CHOOSE_MODEL_PART_NAME",
                "constrained"                    : "fixed",
                "friction_coefficient"           : 0,
                "penalty_coefficient"            : 0
            }  """ )

        # admissible values for "constrained":
        # - "fixed"
        # - "contact"
        # - "slip"
        # - "contact_slip"
        # - "friction"

        settings.ValidateAndAssignDefaults(default_parameters)

        self.model_part = Model[settings["model_part_name"].GetString()]

        # get friction parameters
        self.friction_coefficient = settings["friction_coefficient"].GetDouble()
        self.penalty_coefficient = settings["penalty_coefficient"].GetDouble()

        # Compute the normal on the nodes of interest -
        # Note that the model part employed here is supposed to only have slip "conditions"
        KratosMultiphysics.NormalCalculationUtils().CalculateOnSimplex(self.model_part, self.model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE])

        self.contact = (settings["constrained"].GetString() == "contact")

        for condition in self.model_part.Conditions:
            condition.Set(KratosMultiphysics.SLIP, True)
            condition.Set(KratosMultiphysics.CONTACT, self.contact)
            condition.SetValue(KratosMultiphysics.FRICTION_COEFFICIENT, self.friction_coefficient)
            condition.SetValue(KratosMultiphysics.PENALTY_COEFFICIENT, self.penalty_coefficient)

        for node in self.model_part.Nodes:
            node.Set(KratosMultiphysics.SLIP, True)
            node.SetValue(KratosMultiphysics.FRICTION_COEFFICIENT, self.friction_coefficient)
            node.SetValue(KratosMultiphysics.PENALTY_COEFFICIENT, self.penalty_coefficient)
