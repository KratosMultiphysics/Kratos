from KratosMultiphysics import *
from KratosMultiphysics.SolidMechanicsApplication import *

def Factory(settings, Model):
    if(type(settings) != Parameters):
        raise Exception("Expected input shall be a Parameters object, encapsulating a json string")
    return ApplyScalarOnConditionsProcess(Model, settings["Parameters"])

##all the processes python processes should be derived from "python_process"

class ApplyScalarOnConditionsProcess(Process):
    def __init__(self, Model, settings ):
        Process.__init__(self)        
        
        self.model_part = Model[settings["model_part_name"].GetString()]
        self.var = globals()[settings["variable_name"].GetString()]
        
        #check if variable type is a scalar or a component
        if(type(self.var) != type(DISPLACEMENT_X) and type(self.var) != type(PRESSURE)):
            raise Exception("Variable type is incorrect. Must be a scalar or a component")
        
        self.value = settings["value"].GetDouble()
        
    def ExecuteInitialize(self):
        for cond in self.model_part.Conditions:
            cond.SetValue(self.var,self.value)

    
        
    
