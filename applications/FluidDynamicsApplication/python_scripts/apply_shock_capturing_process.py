import KratosMultiphysics
import KratosMultiphysics.FluidDynamicsApplication as KratosFluid

def Factory(settings, Model):
    if(type(settings) != KratosMultiphysics.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return ApplyShockCapturingProcess(Model, settings["Parameters"])

class ApplyShockCapturingProcess(KratosMultiphysics.Process):

    def __init__(self, Model, settings):

        KratosMultiphysics.Process.__init__(self)

        # Check the default values
        settings = KratosMultiphysics.Parameters("""{
            "model_part_name" : "FluidModelPart",
            "calculate_nodal_area_at_each_step" : false,
            "shock_sensor" : true,
            "shear_sensor" : true,
            "thermal_sensor" : false,
            "thermally_coupled_formulation" : false
        }""")

        settings.ValidateAndAssignDefaults(default_settings)

    def Initialize(self):
        self.shock_process = KratosCFD.ShockCapturingProcess(self.model, settings)
        self.shock_process.Check()
        self.shock_process.ExecuteInitialize()
    
    def FinalizeSolutionStep(self):
        super().FinalizeSolutionStep()
        self.shock_process.ExecuteFinalizeSolutionStep()
