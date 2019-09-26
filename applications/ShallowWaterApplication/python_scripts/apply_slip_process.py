import KratosMultiphysics as KM

def Factory(settings, Model):
    if not isinstance(settings, KM.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return ApplySlipProcess(Model, settings["Parameters"])

## All the processes python should be derived from "Process"
class ApplySlipProcess(KM.Process):
    def __init__(self, Model, settings ):
        KM.Process.__init__(self)

        default_parameters = KM.Parameters("""
            {
                "model_part_name":"PLEASE_CHOOSE_MODEL_PART_NAME"
            }
            """
            )
        settings.ValidateAndAssignDefaults(default_parameters)

        self.model_part = Model[settings["model_part_name"].GetString()]
        self.domain_size = self.model_part.ProcessInfo[KM.DOMAIN_SIZE]

    def ExecuteInitialize(self):
        # Compute the normal on the nodes of interest
        KM.NormalCalculationUtils().CalculateOnSimplex(self.model_part, self.domain_size)

        #TODO: Remove the IS_STRUCTURE variable set as soon as the flag SLIP migration is done
        for node in self.model_part.Nodes:
            node.Set(KM.SLIP)
            node.SetValue(KM.IS_STRUCTURE,1.0)
            node.SetSolutionStepValue(KM.IS_STRUCTURE,0,1.0)
            node.SetSolutionStepValue(KM.MESH_VELOCITY,0,[0.0, 0.0, 0.0])

    def ExecuteInitializeSolutionStep(self):
        KM.VariableUtils().SetScalarVar(KM.IS_STRUCTURE, 0.0, self.model_part.Nodes)
        KM.VariableUtils().SetFlag(KM.SLIP, False, self.model_part.Nodes)
        KM.VariableUtils().SetVectorVar(KM.MESH_VELOCITY, [0.0, 0.0, 0.0], self.model_part.Nodes)
