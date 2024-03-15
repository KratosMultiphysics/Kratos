import KratosMultiphysics

def Factory(settings: KratosMultiphysics.Parameters, model: KratosMultiphysics.Model) -> KratosMultiphysics.OutputProcess:
    if not isinstance(model, KratosMultiphysics.Model):
        raise Exception("expected input shall be a Model object, encapsulating a json string")
    if not isinstance(settings, KratosMultiphysics.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return ApplyMPMFrictionBoundaryProcess(model, settings["Parameters"])

class ApplyMPMFrictionBoundaryProcess(KratosMultiphysics.Process):

    def __init__(self, model: KratosMultiphysics.Model, settings: KratosMultiphysics.Parameters) -> None:
        super().__init__()

        self.settings = settings
        self.settings.ValidateAndAssignDefaults(self.GetDefaultParameters())

        self.model_part = model[settings["model_part_name"].GetString()]
        domain_size = self.model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE]

        # Note that here we are assuming a unique condition geometry in model part
        pts_num = 0
        for condition in self.model_part.Conditions:
            pts_num = condition.GetGeometry().PointsNumber()
            break

        # Get the condition registry name
        cond_reg_name = f"MPMGridNavierSlipFrictionCondition{domain_size}D{pts_num}N"

        # Substitute current conditions by the corresponding ones implementing the wall model
        replace_settings = KratosMultiphysics.Parameters("""{
            "element_name" : "",
            "condition_name" : ""
        }""")
        replace_settings["condition_name"].SetString(cond_reg_name)
        KratosMultiphysics.ReplaceElementsAndConditionsProcess(self.model_part, replace_settings).Execute()

        KratosMultiphysics.VariableUtils().SetFlag(KratosMultiphysics.SLIP, True, self.model_part.Nodes)
        KratosMultiphysics.VariableUtils().SetNonHistoricalVariable(
            KratosMultiphysics.FRICTION_COEFFICIENT,
            settings["friction_coefficient"].GetDouble(),
            self.model_part.Conditions)

        for node in self.model_part.Nodes:
            # Needed to distinguish between conforming and non-conforming BCs (penalty/lagrange)
            node.SetValue(KratosMultiphysics.IS_STRUCTURE,1.0)
            node.SetSolutionStepValue(KratosMultiphysics.IS_STRUCTURE,0,1.0)

        # Compute the normal on the nodes of interest
        KratosMultiphysics.NormalCalculationUtils().CalculateOnSimplex(self.model_part, domain_size)

    @staticmethod
    def GetDefaultParameters() -> KratosMultiphysics.Parameters:
        """Return the default parameters."""
        return KratosMultiphysics.Parameters("""{
            "model_part_name" : "PLEASE_CHOOSE_MODEL_PART_NAME",
            "calculate_normals_at_each_step" : false,
            "friction_coefficient" : 0.0
        }""")

    def ExecuteInitializeSolutionStep(self) -> None:
        """Recompute the normals if needed"""
        if self.settings["calculate_normals_at_each_step"].GetBool():
            KratosMultiphysics.NormalCalculationUtils().CalculateOnSimplex(self.model_part, self.model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE])
