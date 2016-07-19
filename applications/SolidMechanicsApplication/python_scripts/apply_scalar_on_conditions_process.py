import KratosMultiphysics

def Factory(settings, Model):
    if(type(settings) != KratosMultiphysics.Parameters):
        raise Exception("Expected input shall be a Parameters object, encapsulating a json string")
    return ApplyScalarOnConditionsProcess(Model, settings["Parameters"])

##all the processes python processes should be derived from "python_process"

class ApplyScalarOnConditionsProcess(KratosMultiphysics.Process):
    def __init__(self, Model, settings ):
        KratosMultiphysics.Process.__init__(self)        
        
        self.model_part = Model[settings["model_part_name"].GetString()]
        self.var = KratosMultiphysics.KratosGlobals.GetVariable(settings["variable_name"].GetString())
        
        #check if variable type is a scalar or a component
        if(type(self.var) != KratosMultiphysics.Array1DComponentVariable and type(self.var) != KratosMultiphysics.DoubleVariable):
            raise Exception("Variable type is incorrect. Must be a scalar or a component")
        
        self.value = settings["value"].GetDouble()
        
    def ExecuteInitialize(self):
        for cond in self.model_part.Conditions:
            cond.SetValue(self.var,self.value)

    
        
    
