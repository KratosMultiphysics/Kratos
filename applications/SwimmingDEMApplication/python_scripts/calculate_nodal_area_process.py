import KratosMultiphysics

def Factory(settings, Model):
    if not isinstance(settings, KratosMultiphysics.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")

    if settings["Parameters"]["fixed_mesh"].GetBool():
        return CalculateNodalAreaProcess(Model, settings["Parameters"])
    else:
        return CalculateNodalAreaWithRemeshingProcess(Model, settings["Parameters"])

class CalculateNodalAreaProcess(KratosMultiphysics.Process):
    def __init__(self, Model, settings):
        super().__init__()
        default_settings = KratosMultiphysics.Parameters("""
            {
                "model_part_name" : "FluidModelPart.Parts_fluid",
                "domain_size" : 3,
                "fixed_mesh": false
            }
            """)

        settings.ValidateAndAssignDefaults(default_settings)
        self.fluid_model_part = Model[settings["model_part_name"].GetString()]
        dimension = settings["domain_size"].GetInt()
        self.area_calculator = KratosMultiphysics.CalculateNodalAreaProcess(self.fluid_model_part, dimension)

    def ExecuteInitialize(self):
        self.area_calculator.Execute()

    def ExecuteInitializeSolutionStep(self):
        pass

class CalculateNodalAreaWithRemeshingProcess(CalculateNodalAreaProcess):
    def __init__(self, Model, settings):
        super(CalculateNodalAreaWithRemeshingProcess, self).__init__(Model, settings)

    def ExecuteInitializeSolutionStep(self):
        self.area_calculator.Execute()