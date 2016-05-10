from KratosMultiphysics import *
from KratosMultiphysics.SolidMechanicsApplication import *

def Factory(settings, Model):
    if(type(settings) != Parameters):
        raise Exception("Expected input shall be a Parameters object, encapsulating a json string")
    return ApplyVectorOnConditionsProcess(Model, settings["Parameters"])

##all the processes python processes should be derived from "python_process"

class ApplyVectorOnConditionsProcess(Process):
    def __init__(self, Model, settings ):
        Process.__init__(self)
        
        self.model_part = Model[settings["model_part_name"].GetString()]
        self.var = globals()[settings["variable_name"].GetString()]

        #check if variable type is a vector
        if type(self.var) != type(DISPLACEMENT):
            raise Exception("Variable type is incorrect. Must be a three-component vector.")
        
        self.value = [settings["value"][0].GetDouble(),settings["value"][1].GetDouble(),settings["value"][2].GetDouble()]
        
    def ExecuteInitialize(self):
        for cond in self.model_part.Conditions:
            cond.SetValue(self.var,self.value)

    
        
    
