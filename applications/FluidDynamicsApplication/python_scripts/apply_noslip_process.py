from KratosMultiphysics import * 

def Factory(settings, Model):
    if(type(settings) != Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return ApplyNoSlipProcess(Model, settings["Parameters"])

##all the processes python processes should be derived from "python_process"
class ApplyNoSlipProcess(Process):
    def __init__(self, Model, settings):
        
        Process.__init__(self)
        
        model_part = Model[settings["model_part_name"].GetString()]
        settings.AddEmptyValue("is_fixed").SetBool(True)
        settings.AddEmptyValue("value").SetDouble(0)
        
        # For X component
        settings.AddEmptyValue("variable_name").SetString("VELOCITY_X")
        self.x = ApplyConstantScalarValueProcess(model_part,settings)
        
        # For Y component
        settings.AddEmptyValue("variable_name").SetString("VELOCITY_Y")
        self.y = ApplyConstantScalarValueProcess(model_part,settings)
        
        #For Z component
        settings.AddEmptyValue("variable_name").SetString("VELOCITY_Z")
        self.z = ApplyConstantScalarValueProcess(model_part,settings)

    def ExecuteInitialize(self):
        
        self.x.ExecuteInitialize()
        self.y.ExecuteInitialize()
        self.z.ExecuteInitialize()
