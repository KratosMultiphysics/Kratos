import KratosMultiphysics as KM

def Factory(settings, Model):
    if not isinstance(settings, KM.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return ApplySlipProcess(Model, settings["Parameters"])

## All the processes python should be derived from "Process"
class ApplySlipProcess(KM.Process):
    def __init__(self, Model, settings ):
        KM.Process.__init__(self)

        default_parameters = KM.Parameters("""{
                "model_part_name" : "PLEASE_CHOOSE_MODEL_PART_NAME"
            }""")
        settings.ValidateAndAssignDefaults(default_parameters)

        self.model_part = Model[settings["model_part_name"].GetString()]
        self.domain_size = self.model_part.ProcessInfo[KM.DOMAIN_SIZE]

    def ExecuteInitialize(self):
        KM.NormalCalculationUtils().CalculateOnSimplex(self.model_part, self.domain_size)
        KM.VariableUtils().SetFlag(KM.SLIP, False, self.model_part.Nodes)
        KM.VariableUtils().SetVectorVar(KM.MESH_VELOCITY, [0.0, 0.0, 0.0], self.model_part.Nodes)

    def ExecuteInitializeSolutionStep(self):
        KM.VariableUtils().SetFlag(KM.SLIP, False, self.model_part.Nodes)
        KM.VariableUtils().SetVectorVar(KM.MESH_VELOCITY, [0.0, 0.0, 0.0], self.model_part.Nodes)
