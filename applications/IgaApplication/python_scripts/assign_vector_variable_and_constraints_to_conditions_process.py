import KratosMultiphysics

def Factory(settings, model):
    if not isinstance(settings, KratosMultiphysics.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return AssignVectorVariableAndConstraintsToConditionProcess(model, settings["Parameters"])

## All the processes python should be derived from "Process"
class AssignVectorVariableAndConstraintsToConditionProcess(KratosMultiphysics.Process):
    def __init__(self, model, settings ):
        KratosMultiphysics.Process.__init__(self)

        default_settings = KratosMultiphysics.Parameters("""{
            "model_part_name"      : "please_specify_model_part_name",
            "variable_name"        : "SPECIFY_VARIABLE_NAME",
            "interval"             : [0.0, 1e30],
            "value"                : [10.0, "3*t", "x+y"],
            "local_axes"           : {}
        }""")

        # detect "End" as a tag and replace it by a large number
        if(settings.Has("interval")):
            if(settings["interval"][1].IsString() ):
                if(settings["interval"][1].GetString() == "End"):
                    settings["interval"][1].SetDouble(1e30)
                else:
                    raise Exception("the second value of interval can be \"End\" or a number, interval currently:"+settings["interval"].PrettyPrintJsonString())

        settings.ValidateAndAssignDefaults(default_settings)

        variable_name = settings["variable_name"].GetString()
        variable_type = KratosMultiphysics.KratosGlobals.GetVariableType(variable_name)
        self.variable = KratosMultiphysics.KratosGlobals.GetVariable(variable_name)

        if not variable_type == "Array" and not variable_type == "Vector":
            err_msg  = 'Variable-type of variable: "' + variable_name + '" is incorrect ("'
            err_msg += variable_type + '"). Must be a "Array" or "Vector"'
            raise Exception(msg)

        self.model_part = model[settings["model_part_name"].GetString()]

        self.aux_processes = []

        from KratosMultiphysics.assign_scalar_variable_to_conditions_process import AssignScalarVariableToConditionsProcess

        # loop over components X, Y and Z
        for index, variable_dir in enumerate(["_X", "_Y", "_Z"]):
            variable_component_name = variable_name + variable_dir
            flag_name = "KratosMultiphysics.IgaApplication.IgaFlags.FIX_" + variable_component_name
            if not settings["value"][index].IsNull():
                component_params = KratosMultiphysics.Parameters("{}")
                component_params.AddValue("model_part_name", settings["model_part_name"])
                component_params.AddValue("interval", settings["interval"])
                component_params.AddValue("value",settings["value"][index])
                component_params.AddEmptyValue("variable_name").SetString(variable_name + variable_dir)
                component_params.AddValue("local_axes",settings["local_axes"])
                self.aux_processes.append( AssignScalarVariableToConditionsProcess(model, component_params) )
                KratosMultiphysics.VariableUtils().SetFlag(eval(flag_name), True, self.model_part.Conditions)
            else:
                KratosMultiphysics.VariableUtils().SetFlag(eval(flag_name), False, self.model_part.Conditions)

    def ExecuteInitialize(self):
        for process in self.aux_processes:
            process.ExecuteInitialize()

    def ExecuteBeforeSolutionLoop(self):
        for process in self.aux_processes:
            process.ExecuteBeforeSolutionLoop()

    def ExecuteInitializeSolutionStep(self):
        for process in self.aux_processes:
            process.ExecuteInitializeSolutionStep()

    def ExecuteBeforeOutputStep(self):
        for process in self.aux_processes:
            process.ExecuteBeforeOutputStep()

    def ExecuteAfterOutputStep(self):
        for process in self.aux_processes:
            process.ExecuteAfterOutputStep()

    def ExecuteFinalizeSolutionStep(self):
        for process in self.aux_processes:
            process.ExecuteFinalizeSolutionStep()

    def ExecuteFinalize(self):
        for process in self.aux_processes:
            process.ExecuteFinalize()

    def Check(self):
        for process in self.aux_processes:
            process.Check()
