import KratosMultiphysics
from math import *

def Factory(settings, current_model):
    if(type(settings) != KratosMultiphysics.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return AssignVectorVariableProcess(current_model, settings["Parameters"])

## All the processes python should be derived from "Process"
class AssignVectorVariableProcess(KratosMultiphysics.Process):
    def __init__(self, current_model, settings ):
        KratosMultiphysics.Process.__init__(self)

        default_settings = KratosMultiphysics.Parameters("""
            {
                "help"            : "This process assigns a given value (vector) to the nodes belonging a certain submodelpart",
                "mesh_id"              : 0,
                "model_part_name"      : "please_specify_model_part_name",
                "variable_name"        : "SPECIFY_VARIABLE_NAME",
                "interval"             : [0.0, 1e30],
                "value"                : [10.0, "3*t", "x+y"],
                "constrained"          : [true,true,true],
                "local_axes"           : {}
            }
            """
            )
        #example of admissible values for "value" : [10.0, "3*t", "x+y"]

        ##trick to ensure that if someone sets constrained as a single bool, it is transformed to a vector
        if(settings.Has("constrained")):
            if(settings["constrained"].IsBool()):
                is_fixed = settings["constrained"].GetBool()
                #print("is_fixed = ",is_fixed)
                settings["constrained"] = default_settings["constrained"]
                for i in range(3):
                    settings["constrained"][i].SetBool(is_fixed)

        settings.ValidateAndAssignDefaults(default_settings)

        self.variable = KratosMultiphysics.KratosGlobals.GetVariable(settings["variable_name"].GetString())
        if(type(self.variable) != KratosMultiphysics.Array1DVariable3 and type(self.variable) != KratosMultiphysics.VectorVariable):
            msg = "Error in AssignVectorVariableProcess. Variable type of variable : " + settings["variable_name"].GetString() + " is incorrect . Must be a vector or array3"
            raise Exception(msg)

        self.model_part = current_model[settings["model_part_name"].GetString()]

        self.aux_processes = []

        import assign_scalar_variable_process

        # loop over components X, Y and Z
        for indice,variable in enumerate(["_X", "_Y", "_Z"]):
            if(not settings["value"][indice].IsNull()):
                i_params = KratosMultiphysics.Parameters("{}")
                i_params.AddValue("model_part_name",settings["model_part_name"])
                i_params.AddValue("mesh_id",settings["mesh_id"])
                i_params.AddEmptyValue("constrained").SetBool(settings["constrained"][indice].GetBool())
                i_params.AddValue("interval",settings["interval"])
                i_params.AddValue("value",settings["value"][indice])
                i_params.AddEmptyValue("variable_name").SetString(settings["variable_name"].GetString() + variable)
                i_params.AddValue("local_axes",settings["local_axes"])
                self.aux_processes.append( assign_scalar_variable_process.AssignScalarVariableProcess(Model, i_params) )

    def ExecuteBeforeSolutionLoop(self):
        self.ExecuteInitializeSolutionStep()

    def ExecuteInitializeSolutionStep(self):
        for process in self.aux_processes:
            process.ExecuteInitializeSolutionStep()

    def ExecuteFinalizeSolutionStep(self):
        for process in self.aux_processes:
            process.ExecuteFinalizeSolutionStep()
