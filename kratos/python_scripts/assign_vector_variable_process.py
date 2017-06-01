import KratosMultiphysics
from math import *

def Factory(settings, Model):
    if(type(settings) != KratosMultiphysics.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return AssignVectorVariableProcess(Model, settings["Parameters"])

##all the processes python processes should be derived from "python_process"
class AssignVectorVariableProcess(KratosMultiphysics.Process):
    def __init__(self, Model, settings ):
        KratosMultiphysics.Process.__init__(self)

        default_settings = KratosMultiphysics.Parameters("""
            {
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


        #detect "End" as a tag and replace it by a large number
        if(settings.Has("interval")):
            if(settings["interval"][1].IsString() ):
                if(settings["interval"][1].GetString() == "End"):
                    settings["interval"][1].SetDouble(1e30) # = default_settings["interval"][1]
                else:
                    raise Exception("the second value of interval can be \"End\" or a number, interval currently:"+settings["interval"].PrettyPrintJsonString())

        #print(settings.PrettyPrintJsonString())

        settings.ValidateAndAssignDefaults(default_settings)
        
        self.variable = KratosMultiphysics.KratosGlobals.GetVariable(settings["variable_name"].GetString())
        if(type(self.variable) != KratosMultiphysics.Array1DVariable3 and type(self.variable) != KratosMultiphysics.VectorVariable):
            msg = "Error in AssignVectorVariableProcess. Variable type of variable : " + settings["variable_name"].GetString() + " is incorrect . Must be a vector or array3"
            raise Exception(msg)

        self.model_part = Model[settings["model_part_name"].GetString()]

        self.aux_processes = []

        import assign_scalar_variable_process

        #component X
        if(not settings["value"][0].IsNull()):
            x_params = KratosMultiphysics.Parameters("{}")
            x_params.AddValue("model_part_name",settings["model_part_name"])
            x_params.AddValue("mesh_id",settings["mesh_id"])
            (x_params.AddEmptyValue("constrained")).SetBool(settings["constrained"][0].GetBool())
            x_params.AddValue("interval",settings["interval"])
            x_params.AddValue("value",settings["value"][0])
            x_params.AddEmptyValue("variable_name").SetString(settings["variable_name"].GetString() + "_X")
            x_params.AddValue("local_axes",settings["local_axes"])
            self.aux_processes.append( assign_scalar_variable_process.AssignScalarVariableProcess(Model, x_params) )

        #component Y
        if(not settings["value"][1].IsNull()):
            y_params = KratosMultiphysics.Parameters("{}")
            y_params.AddValue("model_part_name",settings["model_part_name"])
            y_params.AddValue("mesh_id",settings["mesh_id"])
            y_params.AddEmptyValue("constrained").SetBool(settings["constrained"][1].GetBool())
            y_params.AddValue("interval",settings["interval"])
            y_params.AddValue("value",settings["value"][1])
            y_params.AddEmptyValue("variable_name").SetString(settings["variable_name"].GetString() + "_Y")
            y_params.AddValue("local_axes",settings["local_axes"])
            self.aux_processes.append( assign_scalar_variable_process.AssignScalarVariableProcess(Model, y_params) )

        #component Z
        if(not settings["value"][2].IsNull()):
            z_params = KratosMultiphysics.Parameters("{}")
            z_params.AddValue("model_part_name",settings["model_part_name"])
            z_params.AddValue("mesh_id",settings["mesh_id"])
            z_params.AddEmptyValue("constrained").SetBool(settings["constrained"][2].GetBool())
            z_params.AddValue("interval",settings["interval"])
            z_params.AddValue("value",settings["value"][2])
            z_params.AddEmptyValue("variable_name").SetString(settings["variable_name"].GetString() + "_Z")
            z_params.AddValue("local_axes",settings["local_axes"])
            self.aux_processes.append( assign_scalar_variable_process.AssignScalarVariableProcess(Model, z_params) )

        # print("Finished construction of AssignVectorProcess Process")

    def ExecuteInitializeSolutionStep(self):
        for process in self.aux_processes:
            process.ExecuteInitializeSolutionStep()

    def ExecuteFinalizeSolutionStep(self):
        #print("---")
        for process in self.aux_processes:
            #print("current_time = ", self.model_part.ProcessInfo[KratosMultiphysics.TIME], " interval = ", process.interval)
            process.ExecuteFinalizeSolutionStep()
