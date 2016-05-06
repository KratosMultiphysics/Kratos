from KratosMultiphysics import *

## This proces sets the value of a scalar variable using the ApplyConstantScalarValueProcess.
## In this case, the scalar value is automatically fixed.

def Factory(settings, Model):
    if(type(settings) != Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return ImposeScalarValueProcess(Model, settings["Parameters"])

## All the processes python processes should be derived from "python_process"
class ImposeScalarValueProcess(ApplyConstantVectorValueProcess):
    def __init__(self, Model, settings ):
        Process.__init__(self)
        
        model_part = Model[settings["model_part_name"].GetString()]  
        
        settings.AddEmptyValue("is_fixed").SetBool(True)       
        
        ApplyConstantScalarValueProcess.__init__(self,model_part, settings)


        

            

        
        
        
