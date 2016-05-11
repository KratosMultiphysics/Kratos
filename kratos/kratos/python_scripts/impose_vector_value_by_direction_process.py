from KratosMultiphysics import *

## This proces sets the value of a vector variable using the ApplyConstantVectorValueProcess.
## In this case, the three components of the vector are automatically fixed.

def Factory(settings, Model):
    if(type(settings) != Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return ImposeVectorValueByDirectionProcess(Model, settings["Parameters"])

## All the processes python processes should be derived from "python_process"
class ImposeVectorValueByDirectionProcess(ApplyConstantVectorValueProcess):
    def __init__(self, Model, settings ):
        Process.__init__(self)
        
        model_part = Model[Parameters["model_part_name"].GetString()]  
        
        settings.AddValue("is_fixed_x",True)        
        settings.AddValue("is_fixed_y",True)        
        settings.AddValue("is_fixed_z",True)
        
        ApplyConstantVectorValueProcess.__init__(self,model_part, settings)


        

            

        
        
        
