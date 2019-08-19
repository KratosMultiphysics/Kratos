import KratosMultiphysics as KratosMultiphysics
import KratosMultiphysics.StructuralMechanicsApplication as StructuralMechanicsApplication

def Factory(settings, Model):
    if(type(settings) != KratosMultiphysics.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return ApplyWeakSlidingProcess(Model, settings["Parameters"])



class ApplyWeakSlidingProcessStruct(KratosMultiphysics.Process):

    def __init__(self, Model, settings ):
        KratosMultiphysics.Process.__init__(self)
        default_settings = KratosMultiphysics.Parameters("""
        {
            "model_part_name_slave"           : "example_part_slave",
            "model_part_name_master"          : "example_part_master",
            "computing_model_part_name"       : "computing_domain",
            "element_id"                      : 1,
            "property_id"                     : 1,
            "debug_info"                      : false
        }
        """)
        default_settings.ValidateAndAssignDefaults(settings)

        # The computing model part
        self.computing_model_part = Model["Structure"].GetSubModelPart(settings["computing_model_part_name"].GetString())
        self.wip_model_part      = Model["Structure"].GetSubModelPart(settings["model_part_name_slave"].GetString())

        self.custom_process = StructuralMechanicsApplication.ApplyWeakSlidingProcessStruct(Model["Structure"], settings)


    def ExecuteInitialize(self):
        self.custom_process.ExecuteInitialize()

    def ExecuteInitializeSolutionStep(self):
        self.custom_process.ExecuteInitializeSolutionStep()

    def ExecuteFinalizeSolutionStep(self):
        self.custom_process.ExecuteFinalizeSolutionStep()