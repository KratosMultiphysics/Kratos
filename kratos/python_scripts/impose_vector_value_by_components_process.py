from KratosMultiphysics import *

## This proces sets the value of a vector variable component-by-component.
## In this case, the fixicity is given by the user and some of the components may not be fixed.

def Factory(settings, Model):
    if(type(settings) != Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return ImposeVectorValueByComponentsProcess(Model, settings["Parameters"])

## All the processes python processes should be derived from "python_process"
class ImposeVectorValueByComponentsProcess(Process):
    def __init__(self, Model, settings ):
        Process.__init__(self)
        
        model_part = Model[settings["model_part_name"].GetString()]
        variable_name = settings["variable_name"].GetString()
        
        # Auxiliar vector with the fixicity settings    
        self.fixicity_vec = [settings["is_fixed_x"].GetBool(),
                             settings["is_fixed_y"].GetBool(),
                             settings["is_fixed_z"].GetBool()]

        if settings["is_fixed_x"].GetBool() == True:
            # Auxiliar x-component parameters creation
            x_params = Parameters("{}")           
            x_params.AddValue("model_part_name",settings["model_part_name"])
            x_params.AddValue("mesh_id",settings["mesh_id"])
            x_params.AddValue("is_fixed",settings["is_fixed_x"])
            x_params.AddValue("value",settings["value"][0])
            x_params.AddEmptyValue("variable_name").SetString(variable_name+"_X")
            
            self.x_component_process = ApplyConstantScalarValueProcess(model_part, x_params)
            
        if settings["is_fixed_y"].GetBool() == True:
            # Auxiliar y-component parameters creation
            y_params = Parameters("{}")
            y_params.AddValue("model_part_name",settings["model_part_name"])
            y_params.AddValue("mesh_id",settings["mesh_id"])
            y_params.AddValue("is_fixed",settings["is_fixed_y"])
            y_params.AddValue("value",settings["value"][1])
            y_params.AddEmptyValue("variable_name").SetString(variable_name+"_Y")
            
            self.y_component_process = ApplyConstantScalarValueProcess(model_part, y_params)
            
        if settings["is_fixed_z"].GetBool() == True:
            # Auxiliar x-component parameters creation
            z_params = Parameters("{}")
            z_params.AddValue("model_part_name",settings["model_part_name"])
            z_params.AddValue("mesh_id",settings["mesh_id"])
            z_params.AddValue("is_fixed",settings["is_fixed_z"])
            z_params.AddValue("value",settings["value"][2])
            z_params.AddEmptyValue("variable_name").SetString(variable_name+"_Z")
            
            self.z_component_process = ApplyConstantScalarValueProcess(model_part, z_params)
                    
        
    def ExecuteInitialize(self):
               
        if self.fixicity_vec[0] == True:
            self.x_component_process.ExecuteInitialize()
            
        if self.fixicity_vec[1] == True:
            self.y_component_process.ExecuteInitialize()
            
        if self.fixicity_vec[2] == True:
            self.z_component_process.ExecuteInitialize()
        
        

            

        
        
        
