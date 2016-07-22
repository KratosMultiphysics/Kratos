from KratosMultiphysics import * 

def Factory(settings, Model):
    if(type(settings) != Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return SetInterfaceProcess(Model, settings["Parameters"])

##all the processes python processes should be derived from "python_process"
class SetInterfaceProcess(Process):
    def __init__(self, Model, settings):
        
        Process.__init__(self)
        
        model_part = Model[settings["model_part_name"].GetString()]
        settings.AddEmptyValue("value").SetDouble(1.0)
        settings.AddEmptyValue("is_fixed").SetBool(True)
        
        self.set_interface_process = ApplyConstantScalarValueProcess(model_part,settings)
        
        
    def ExecuteInitialize(self):
        
        self.set_interface_process.ExecuteInitialize()
