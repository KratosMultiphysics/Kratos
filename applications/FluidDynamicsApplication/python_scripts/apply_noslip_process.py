import KratosMultiphysics

def Factory(settings, Model):
    if(type(settings) != KratosMultiphysics.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return ApplyNoSlipProcess(Model, settings["Parameters"])


class ApplyNoSlipProcess(KratosMultiphysics.Process):
    def __init__(self, Model, settings):

        KratosMultiphysics.Process.__init__(self)

        default_settings = KratosMultiphysics.Parameters("""
            {
                "mesh_id"              : 0,
                "model_part_name"      : "please_specify_model_part_name",
                "interval"             : [0.0, "End"]
            }
            """
            )

        # Detect "End" as a tag and replace it by a large number
        if(settings.Has("interval")):
            if(settings["interval"][1].IsString() ):
                if(settings["interval"][1].GetString() == "End"):
                    settings["interval"][1].SetDouble(1e30) # = default_settings["interval"][1]
                else:
                    raise Exception("the second value of interval can be \"End\" or a number, interval currently:"+settings["interval"].PrettyPrintJsonString())

        settings.ValidateAndAssignDefaults(default_settings)

        self.aux_processes = []

        import experimental_assign_value_process

        # Component X
        x_params = KratosMultiphysics.Parameters("{}")
        x_params.AddValue("model_part_name",settings["model_part_name"])
        x_params.AddValue("mesh_id",settings["mesh_id"])
        x_params.AddEmptyValue("constrained").SetBool(True)
        x_params.AddValue("interval",settings["interval"])
        x_params.AddEmptyValue("value").SetDouble(0.0)
        x_params.AddEmptyValue("variable_name").SetString("VELOCITY_X")
        self.aux_processes.append( experimental_assign_value_process.AssignValueProcess(Model, x_params) )

        # Component Y
        y_params = KratosMultiphysics.Parameters("{}")
        y_params.AddValue("model_part_name",settings["model_part_name"])
        y_params.AddValue("mesh_id",settings["mesh_id"])
        y_params.AddEmptyValue("constrained").SetBool(True)
        y_params.AddValue("interval",settings["interval"])
        y_params.AddEmptyValue("value").SetDouble(0.0)
        y_params.AddEmptyValue("variable_name").SetString("VELOCITY_Y")
        self.aux_processes.append( experimental_assign_value_process.AssignValueProcess(Model, y_params) )

        # Component Y
        z_params = KratosMultiphysics.Parameters("{}")
        z_params.AddValue("model_part_name",settings["model_part_name"])
        z_params.AddValue("mesh_id",settings["mesh_id"])
        z_params.AddEmptyValue("constrained").SetBool(True)
        z_params.AddValue("interval",settings["interval"])
        z_params.AddEmptyValue("value").SetDouble(0.0)
        z_params.AddEmptyValue("variable_name").SetString("VELOCITY_Z")
        self.aux_processes.append( experimental_assign_value_process.AssignValueProcess(Model, z_params) )


    def ExecuteInitializeSolutionStep(self):
        for process in self.aux_processes:
            process.ExecuteInitializeSolutionStep()


    def ExecuteFinalizeSolutionStep(self):
        for process in self.aux_processes:
            process.ExecuteFinalizeSolutionStep()
