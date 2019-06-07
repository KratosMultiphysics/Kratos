import KratosMultiphysics
import KratosMultiphysics.RANSModellingApplication as KratosRANS


def Factory(settings, Model):
    if (type(settings) != KratosMultiphysics.Parameters):
        raise Exception(
            "expected input shall be a Parameters object, encapsulating a json string"
        )
    if (type(Model) != KratosMultiphysics.Model):
        raise Exception("expected input shall be a Model object")

    return ApplyCustomBoundaryProcess(Model, settings["Parameters"])


## All the processes python should be derived from "Process"
class ApplyCustomBoundaryProcess(KratosMultiphysics.Process):
    def __init__(self, Model, settings):
        KratosMultiphysics.Process.__init__(self)

        default_parameters = KratosMultiphysics.Parameters("""
            {
                "process_name"              : "PLEASE_SPECIFY_PROCESS_NAME",
                "model_part_name"           : "PLEASE_SPECIFY_MODEL_PART_NAME",
                "is_initialization_process" : false,
                "boundary_condition_type"   : "none",
                "is_constrained"            : true,
                "process_settings"          : {}
            }  """)

        settings.ValidateAndAssignDefaults(default_parameters)
        self.settings = settings
        self.model_part = Model[self.settings["model_part_name"].GetString()]
        self.is_initialization_process = self.settings[
            "is_initialization_process"].GetBool()
        self.boundary_condition_type = self.settings[
            "boundary_condition_type"].GetString()
        self.is_initialized = False

        custom_process_list = [
            "k_turbulent_intensity", "epsilon_turbulent_mixing_length"
        ]
        process_name = self.settings["process_name"].GetString()

        if (not (process_name in custom_process_list)):
            msg = "process_name = \"" + process_name + "\" not found. Valid process names are:\n      "
            msg += "\n      ".join(custom_process_list)
            msg += "\n"
            raise Exception(msg)

        if (not self.is_initialization_process):
            from model_part_factory import ApplyFlagsToNodeList
            ApplyFlagsToNodeList(
                self.model_part.Nodes,
                self.settings["boundary_condition_type"].GetString(), True)
        else:
            self.settings["is_constrained"].SetBool(False)

        if (process_name == "epsilon_turbulent_mixing_length"):
            self.process = KratosRANS.RansEpsilonTurbulentMixingLengthEvaluationProcess(
                self.model_part, self.settings["process_settings"],
                self.settings["is_constrained"].GetBool())
        elif (process_name == "k_turbulent_intensity"):
            self.process = KratosRANS.RansKTurbulentIntensityEvaluationProcess(
                self.model_part, self.settings["process_settings"],
                self.settings["is_constrained"].GetBool())

    def Check(self):
        self.process.Check()

    def ExecuteInitializeSolutionStep(self):
        if (self.model_part.ProcessInfo[KratosRANS.
                                        IS_CO_SOLVING_PROCESS_ACTIVE]):
            if (not self.is_initialized and self.is_initialization_process):
                self.process.Execute()
                self.is_initialized = True
            elif (not self.is_initialization_process):
                self.process.Execute()
