import KratosMultiphysics 
from math import *
        
def Factory(settings, Model):
    if(type(settings) != KratosMultiphysics.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return AssignVectorProcess(Model, settings["Parameters"])

    
##all the processes python processes should be derived from "python_process"
class AssignVectorProcess(KratosMultiphysics.Process):
    def __init__(self, Model, settings ):
        KratosMultiphysics.Process.__init__(self) 

        default_settings = KratosMultiphysics.Parameters("""
            {
                "mesh_id"              : 0,
                "model_part_name"      : "please_specify_model_part_name",
                "variable_name"        : "SPECIFY_VARIABLE_NAME",
                "interval"             : [0.0, 1e30],
                "value"                : [10.0, "3*t", "x+y"],
                "local_axes"           : {}
            }
            """
            )
        
        #example of admissible values for "value" : [10.0, "3*t", "x+y"]
        
        settings.ValidateAndAssignDefaults(default_settings)
        
        self.model_part = Model[settings["model_part_name"].GetString()]
        
        self.aux_processes = []
        
        import experimental_assign_value_process
        
        #component X
        if(not settings["value"][0].IsNull()):
            x_params = KratosMultiphysics.Parameters("{}")           
            x_params.AddValue("model_part_name",settings["model_part_name"])
            x_params.AddValue("mesh_id",settings["mesh_id"])
            x_params.AddEmptyValue("constrained").SetBool(True)
            x_params.AddValue("interval",settings["interval"])
            x_params.AddValue("value",settings["value"][0])
            x_params.AddEmptyValue("variable_name").SetString(settings["variable_name"].GetString() + "_X")
            x_params.AddValue("local_axes",settings["local_axes"])
            self.aux_processes.append( experimental_assign_value_process.AssignValueProcess(Model, x_params) )

        #component Y
        if(not settings["value"][1].IsNull()):
            y_params = KratosMultiphysics.Parameters("{}")           
            y_params.AddValue("model_part_name",settings["model_part_name"])
            y_params.AddValue("mesh_id",settings["mesh_id"])
            y_params.AddEmptyValue("constrained").SetBool(True)
            y_params.AddValue("interval",settings["interval"])
            y_params.AddValue("value",settings["value"][1])
            y_params.AddEmptyValue("variable_name").SetString(settings["variable_name"].GetString() + "_Y")
            y_params.AddValue("local_axes",settings["local_axes"])
            self.aux_processes.append( experimental_assign_value_process.AssignValueProcess(Model, y_params) )

        #component Z
        if(not settings["value"][2].IsNull()):
            z_params = KratosMultiphysics.Parameters("{}")           
            z_params.AddValue("model_part_name",settings["model_part_name"])
            z_params.AddValue("mesh_id",settings["mesh_id"])
            z_params.AddEmptyValue("constrained").SetBool(True)
            z_params.AddValue("interval",settings["interval"])
            z_params.AddValue("value",settings["value"][2])
            z_params.AddEmptyValue("variable_name").SetString(settings["variable_name"].GetString() + "_Z")
            z_params.AddValue("local_axes",settings["local_axes"])
            self.aux_processes.append( experimental_assign_value_process.AssignValueProcess(Model, z_params) )
        
        print("finished construction of AssignVectorProcess Process")
        
    def ExecuteInitializeSolutionStep(self):
        for process in self.aux_processes:
            process.ExecuteInitializeSolutionStep()
            
    def ExecuteFinalizeSolutionStep(self):
        for process in self.aux_processes:
            process.ExecuteInitializeSolutionStep()
