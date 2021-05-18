import KratosMultiphysics as KM

def Factory(settings, model):
    if not isinstance(settings, KM.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return ApplySlipProcess(model, settings["Parameters"])

class ApplySlipProcess(KM.Process):
    """ApplySlipProcess

    This process sets the SLIP flag and computes the
    NORMAL variable on the selected model part
    """
    def __init__(self, model, settings):
        """The constructor of the ApplySlipProcess

        Keyword arguments:
        self -- It signifies an instance of a class.
        model -- The model to be used
        settings -- The ProjectParameters used
        """
        KM.Process.__init__(self)

        default_parameters = KM.Parameters("""{
                "model_part_name" : "PLEASE_CHOOSE_MODEL_PART_NAME",
                "recompute_normals" : false
            }""")
        settings.ValidateAndAssignDefaults(default_parameters)

        self.model_part = model[settings["model_part_name"].GetString()]
        self.recompute_normals = settings["recompute_normals"].GetBool()

    def ExecuteInitialize(self):
        domain_size = self.model_part.ProcessInfo[KM.DOMAIN_SIZE]
        KM.NormalCalculationUtils().CalculateOnSimplex(self.model_part, domain_size)
        KM.VariableUtils().SetFlag(KM.SLIP, True, self.model_part.Nodes)

    def ExecuteInitializeSolutionStep(self):
        if self.recompute_normals:
            self.ExecuteInitialize()
