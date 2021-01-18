import KratosMultiphysics as KratosMultiphysics
import KratosMultiphysics.CableNetApplication as CableNetApplication

def Factory(settings, Model):
    if(type(settings) != KratosMultiphysics.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return ApplyWeakSlidingProcess(Model, settings["Parameters"])



class ApplyWeakSlidingProcess(KratosMultiphysics.Process):

    def __init__(self, Model, settings ):
        KratosMultiphysics.Process.__init__(self)
        default_settings = KratosMultiphysics.Parameters("""
        {
            "model_part_name_slave"           : "example_part_slave",
            "model_part_name_master"          : "example_part_master",
            "computing_model_part_name"       : "Structure",
            "element_id"                      : 1,
            "property_id"                     : 1,
            "debug_info"                      : false
        }
        """)
        default_settings.ValidateAndAssignDefaults(settings)
        self.custom_process = CableNetApplication.ApplyWeakSlidingProcess(Model[settings["computing_model_part_name"].GetString()], settings)


    def ExecuteInitialize(self):
        self.custom_process.ExecuteInitialize()

    def ExecuteInitializeSolutionStep(self):
        self.custom_process.ExecuteInitializeSolutionStep()

    def ExecuteFinalizeSolutionStep(self):
        self.custom_process.ExecuteFinalizeSolutionStep()