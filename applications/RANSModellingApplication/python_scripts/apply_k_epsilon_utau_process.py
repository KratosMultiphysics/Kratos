import KratosMultiphysics
import KratosMultiphysics.RANSModellingApplication as KratosRANS

def Factory(settings, Model):
    if(type(settings) != KratosMultiphysics.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    if(type(Model) != KratosMultiphysics.Model):
        raise Exception("expected input shall be a Model object")

    return ApplyKEpsilonUtauProcess(Model, settings["Parameters"])

## All the processes python should be derived from "Process"
class ApplyKEpsilonUtauProcess(KratosMultiphysics.Process):
    def __init__(self, Model, settings ):
        KratosMultiphysics.Process.__init__(self)

        default_parameters = KratosMultiphysics.Parameters( """
            {
                "model_part_name"           : "PLEASE_CHOOSE_MODEL_PART_NAME",
                "y_plus_model"              : {},
                "constants"                 : {},
                "is_initialization_process" : false,
                "boundary_condition_type"   : "inlet"
            }  """ )

        settings.ValidateAndAssignDefaults(default_parameters)
        self.settings = settings

        self.model_part = Model[settings["model_part_name"].GetString()]
        import rans_y_plus_model_factory

        rans_y_plus_model_factory.InitializeModelPartName(settings["y_plus_model"], Model, self.model_part)

        self.y_plus_model = rans_y_plus_model_factory.Factory(self.settings["y_plus_model"], Model)
        self.settings.RemoveValue("y_plus_model")

        if (not self.settings["is_initialization_process"].GetBool()):
            from model_part_factory import ApplyFlagsToNodeList
            ApplyFlagsToNodeList(self.model_part.Nodes, self.settings["boundary_condition_type"].GetString(), True)
        self.settings.RemoveValue("boundary_condition_type")

        self.k_epsilon_utau_process = KratosRANS.RansKEpsilonEvaluationUtauProcess(Model, self.settings, self.y_plus_model)

    def Check(self):
        self.y_plus_model.Check()
        self.k_epsilon_utau_process.Check()

    def ExecuteInitializeSolutionStep(self):
        if (self.model_part.ProcessInfo[KratosRANS.IS_CO_SOLVING_PROCESS_ACTIVE]):
            self.k_epsilon_utau_process.Execute()
