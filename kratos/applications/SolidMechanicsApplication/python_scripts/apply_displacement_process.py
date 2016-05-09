from KratosMultiphysics import *

#~ def Factory(settings, Model):
    #~ if(type(settings) != Parameters):
        #~ raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    #~ return ApplyDisplacementProcess(Model, settings["Parameters"])

##all the processes python processes should be derived from "python_process"
#~ class ApplyDisplacementProcess(ApplyConstantVectorValueProcess):
    #~ def __init__(self, Model, Parameters ):
        #~ model_part = Model[Parameters["model_part_name"].GetString()]      
        #~ Parameters.AddEmptyValue("variable_name").SetString("DISPLACEMENT")
        #~ # Parameters.AddEmptyValue("is_fixed_x").SetBool(True)
        #~ # Parameters.AddEmptyValue("is_fixed_y").SetBool(True)
        #~ # Parameters.AddEmptyValue("is_fixed_z").SetBool(True)
               
        #~ ApplyConstantVectorValueProcess.__init__(self,model_part, Parameters)

def Factory(settings, Model):
    if(type(settings) != Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return ApplyDisplacementProcess(Model, settings["Parameters"])

class ApplyDisplacementProcess(Process):
    def __init__(self, Model, settings ):
        Process.__init__(self)
        
        model_part = Model[settings["model_part_name"].GetString()]

        # Auxiliar x-component parameters creation
        x_params = Parameters("{}")           
        x_params.AddValue("model_part_name",settings["model_part_name"])
        x_params.AddValue("mesh_id",settings["mesh_id"])
        x_params.AddValue("is_fixed",settings["is_fixed_x"])
        x_params.AddValue("value",settings["direction"][0])
        x_params.AddEmptyValue("variable_name").SetString("DISPLACEMENT_X")
        
        self.x_component_process = ApplyConstantScalarValueProcess(model_part, x_params)
        
        self.model_part = model_part
        
        # Auxiliar y-component parameters creation
        y_params = Parameters("{}")
        y_params.AddValue("model_part_name",settings["model_part_name"])
        y_params.AddValue("mesh_id",settings["mesh_id"])
        y_params.AddValue("is_fixed",settings["is_fixed_y"])
        y_params.AddValue("value",settings["direction"][1])
        y_params.AddEmptyValue("variable_name").SetString("DISPLACEMENT_Y")
        
        self.y_component_process = ApplyConstantScalarValueProcess(model_part, y_params)
            
        # Auxiliar x-component parameters creation
        z_params = Parameters("{}")
        z_params.AddValue("model_part_name",settings["model_part_name"])
        z_params.AddValue("mesh_id",settings["mesh_id"])
        z_params.AddValue("is_fixed",settings["is_fixed_z"])
        z_params.AddValue("value",settings["direction"][2])
        z_params.AddEmptyValue("variable_name").SetString("DISPLACEMENT_Z")
        
        self.z_component_process = ApplyConstantScalarValueProcess(model_part, z_params)
        
    def ExecuteInitialize(self):
        
        self.x_component_process.ExecuteInitialize()
        
        self.y_component_process.ExecuteInitialize()
            
        self.z_component_process.ExecuteInitialize()
