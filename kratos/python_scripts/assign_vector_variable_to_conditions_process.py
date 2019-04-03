import KratosMultiphysics
from math import *

def Factory(settings, Model):
    if(type(settings) != KratosMultiphysics.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return AssignVectorVariableToConditionProcess(Model, settings["Parameters"])

## All the processes python should be derived from "Process"
class AssignVectorVariableToConditionProcess(KratosMultiphysics.Process):
    def __init__(self, Model, settings ):
        KratosMultiphysics.Process.__init__(self)

        default_settings = KratosMultiphysics.Parameters("""
            {
                "help"                 : "This process assigns a given value (vector) to the conditions belonging a certain submodelpart",
                "mesh_id"              : 0,
                "model_part_name"      : "please_specify_model_part_name",
                "variable_name"        : "SPECIFY_VARIABLE_NAME",
                "interval"             : [0.0, 1e30],
                "value"                : [10.0, "3*t", "x+y"],
                "local_axes"           : {}
            }
            """
            )

        # detect "End" as a tag and replace it by a large number
        if(settings.Has("interval")):
            if(settings["interval"][1].IsString() ):
                if(settings["interval"][1].GetString() == "End"):
                    settings["interval"][1].SetDouble(1e30) # = default_settings["interval"][1]
                else:
                    raise Exception("the second value of interval can be \"End\" or a number, interval currently:"+settings["interval"].PrettyPrintJsonString())

        settings.ValidateAndAssignDefaults(default_settings)
        
        self.variable = KratosMultiphysics.KratosGlobals.GetVariable(settings["variable_name"].GetString())
        if(type(self.variable) != KratosMultiphysics.Array1DVariable3 and type(self.variable) != KratosMultiphysics.VectorVariable):
            msg = "Error in AssignVectorToConditionProcess. Variable type of variable : " + settings["variable_name"].GetString() + " is incorrect . Must be a vector or array3"
            raise Exception(msg)

        self.model_part = Model[settings["model_part_name"].GetString()]

        self.aux_processes = []

        import assign_scalar_variable_to_conditions_process

        # component X
        if(not settings["value"][0].IsNull()):
            x_params = KratosMultiphysics.Parameters("{}")
            x_params.AddValue("model_part_name",settings["model_part_name"])
            x_params.AddValue("mesh_id",settings["mesh_id"])
            x_params.AddValue("interval",settings["interval"])
            x_params.AddValue("value",settings["value"][0])
            x_params.AddEmptyValue("variable_name").SetString(settings["variable_name"].GetString() + "_X")
            x_params.AddValue("local_axes",settings["local_axes"])
            self.aux_processes.append( assign_scalar_variable_to_conditions_process.AssignScalarVariableToConditionsProcess(Model, x_params) )

        # component Y
        if(not settings["value"][1].IsNull()):
            y_params = KratosMultiphysics.Parameters("{}")
            y_params.AddValue("model_part_name",settings["model_part_name"])
            y_params.AddValue("mesh_id",settings["mesh_id"])
            y_params.AddValue("interval",settings["interval"])
            y_params.AddValue("value",settings["value"][1])
            y_params.AddEmptyValue("variable_name").SetString(settings["variable_name"].GetString() + "_Y")
            y_params.AddValue("local_axes",settings["local_axes"])
            self.aux_processes.append( assign_scalar_variable_to_conditions_process.AssignScalarVariableToConditionsProcess(Model, y_params) )

        # component Z
        if(not settings["value"][2].IsNull()):
            z_params = KratosMultiphysics.Parameters("{}")
            z_params.AddValue("model_part_name",settings["model_part_name"])
            z_params.AddValue("mesh_id",settings["mesh_id"])
            z_params.AddValue("interval",settings["interval"])
            z_params.AddValue("value",settings["value"][2])
            z_params.AddEmptyValue("variable_name").SetString(settings["variable_name"].GetString() + "_Z")
            z_params.AddValue("local_axes",settings["local_axes"])
            self.aux_processes.append( assign_scalar_variable_to_conditions_process.AssignScalarVariableToConditionsProcess(Model, z_params) )

    def ExecuteInitializeSolutionStep(self):
        for process in self.aux_processes:
            process.ExecuteInitializeSolutionStep()

    def ExecuteFinalizeSolutionStep(self):
        for process in self.aux_processes:
            process.ExecuteFinalizeSolutionStep()
