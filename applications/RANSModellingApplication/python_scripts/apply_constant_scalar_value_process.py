import KratosMultiphysics
import KratosMultiphysics.RANSModellingApplication as KratosRANS

def Factory(settings, Model):
    if(type(settings) != KratosMultiphysics.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return ApplyConstantScalarValueProcess(Model, settings["Parameters"])

## All the processes python should be derived from "Process"
class ApplyConstantScalarValueProcess(KratosMultiphysics.Process):
    def __init__(self, Model, settings ):
        KratosMultiphysics.Process.__init__(self)

        default_parameters = KratosMultiphysics.Parameters( """
            {
                "model_part_name":"PLEASE_CHOOSE_MODEL_PART_NAME",
                "mesh_id": 0,
                "variable_name": "PLEASE_PRESCRIBE_VARIABLE_NAME",
                "is_fixed": false,
                "value" : 1.0,
                "boundary_condition_type" : "structure"
            }  """ )

        settings.ValidateAndAssignDefaults(default_parameters)
        self.settings = settings

        self.model_part = Model[settings["model_part_name"].GetString()]

        from model_part_factory import ApplyFlagsToModelPart
        ApplyFlagsToModelPart(self.model_part, self.settings["boundary_condition_type"].GetString(), True)
        self.settings.RemoveValue("boundary_condition_type")

        self.scalar_value_process = KratosMultiphysics.ApplyConstantScalarValueProcess(self.model_part, self.settings)
        self.is_constant_scalars_applied = False

    def ExecuteInitializeSolutionStep(self):
        if (self.model_part.ProcessInfo[KratosRANS.IS_CO_SOLVING_PROCESS_ACTIVE] and not self.is_constant_scalars_applied):
            self.scalar_value_process.ExecuteInitialize()
            self.is_constant_scalars_applied = True
