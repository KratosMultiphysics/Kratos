import KratosMultiphysics

## This proces sets the value of a vector variable using the ApplyConstantVectorValueProcess.
## In this case, the three components of the vector are automatically fixed.

def Factory(settings, Model):
    if(type(settings) != KratosMultiphysics.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return ImposeVectorValueByDirectionProcess(Model, settings["Parameters"])

## All the processes python processes should be derived from "python_process"
class ImposeVectorValueByDirectionProcess(KratosMultiphysics.ApplyConstantVectorValueProcess):
    def __init__(self, Model, settings ):
        KratosMultiphysics.Process.__init__(self)
        
        model_part = Model[settings["model_part_name"].GetString()]  
        
        settings.AddEmptyValue("is_fixed_x").SetBool(True)
        settings.AddEmptyValue("is_fixed_y").SetBool(True)
        settings.AddEmptyValue("is_fixed_z").SetBool(True)
        
        KratosMultiphysics.ApplyConstantVectorValueProcess.__init__(self,model_part, settings)


        

            

        
        
        
