import KratosMultiphysics as Kratos

def Factory(settings: Kratos.Parameters, model: Kratos.Model):
    if not isinstance(settings, Kratos.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return CheckPointProcess(model, settings["Parameters"])

class CheckPointProcess(Kratos.Process):
    def __init__(self, model: Kratos.Model, settings: Kratos.Parameters):
        Kratos.Process.__init__(self)

        default_settings = Kratos.Parameters("""
        {
            "help"            : "Saves or loads a model part checkpoint using a file serializer",
            "model_part_name" : "please_specify_model_part_name",
            "checkpoint_state": "",
            "file_name"       : ""
        }
        """
        )

        settings.ValidateAndAssignDefaults(default_settings)

        self.model_part = model[settings["model_part_name"].GetString()]
        self.check_point_state = settings["checkpoint_state"].GetString()
        self.file_name = settings["file_name"].GetString()

        if self.check_point_state not in ["save", "load"]:
            raise RuntimeError(f"Unsupported checkpoint state = \"{self.check_point_state}\".Followings are supported:\n\tload\n\tsave")

    def ExecuteInitializeSolutionStep(self):
        if self.check_point_state == "load":
            current_filename = self.file_name.replace("<time>", str(self.model_part.ProcessInfo[Kratos.TIME]))
            current_filename = current_filename.replace("<step>", str(self.model_part.ProcessInfo[Kratos.STEP]))
            Kratos.FileSerializer(current_filename, Kratos.SerializerTraceType.SERIALIZER_NO_TRACE, True).Load("check_point", self.model_part)

    def ExecuteFinalizeSolutionStep(self):
        if self.check_point_state == "save":
            current_filename = self.file_name.replace("<time>", str(self.model_part.ProcessInfo[Kratos.TIME]))
            current_filename = current_filename.replace("<step>", str(self.model_part.ProcessInfo[Kratos.STEP]))
            Kratos.FileSerializer(current_filename, Kratos.SerializerTraceType.SERIALIZER_NO_TRACE, True).Save("check_point", self.model_part)
