from KratosMultiphysics import *
from FluidDynamicsApplication import *

def Factory(settings, Model):
    if(type(settings) != Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return ApplyConstantInletProcess(Model, settings["Parameters"])

##all the processes python processes should be derived from "python_process"
#This process is very similar to apply_custom_velocity_constraints, but fixing every component
class ApplyConstantInletProcess(Process):
    def __init__(self, Model, settings ):
        Process.__init__(self)
        
        model_part = Model[settings["model_part_name"].GetString()]

        ## X DOF
        # Auxiliar x-component parameters creation
        x_params = Parameters("{}")           
        x_params.AddValue("model_part_name",settings["model_part_name"])
        x_params.AddValue("mesh_id",settings["mesh_id"])
        x_params.AddEmptyValue("is_fixed").SetBool(True)
        x_params.AddValue("value",settings["direction"][0])
        x_params.AddEmptyValue("variable_name").SetString("VELOCITY_X")
        
        self.x_component_process = ApplyConstantScalarValueProcess(model_part, x_params)
            
        ## Y DOF
        # Auxiliar y-component parameters creation
        y_params = Parameters("{}")
        y_params.AddValue("model_part_name",settings["model_part_name"])
        y_params.AddValue("mesh_id",settings["mesh_id"])
        y_params.AddEmptyValue("is_fixed").SetBool(True)
        y_params.AddValue("value",settings["direction"][1])
        y_params.AddEmptyValue("variable_name").SetString("VELOCITY_Y")
        
        self.y_component_process = ApplyConstantScalarValueProcess(model_part, y_params)
            
        ## Z DOF
        # Auxiliar x-component parameters creation
        z_params = Parameters("{}")
        z_params.AddValue("model_part_name",settings["model_part_name"])
        z_params.AddValue("mesh_id",settings["mesh_id"])
        z_params.AddEmptyValue("is_fixed").SetBool(True)
        z_params.AddValue("value",settings["direction"][2])
        z_params.AddEmptyValue("variable_name").SetString("VELOCITY_Z")
        
        self.z_component_process = ApplyConstantScalarValueProcess(model_part, z_params)
            
        
    def ExecuteInitialize(self):

        self.x_component_process.ExecuteInitialize()
        self.y_component_process.ExecuteInitialize()
        self.z_component_process.ExecuteInitialize()
    
