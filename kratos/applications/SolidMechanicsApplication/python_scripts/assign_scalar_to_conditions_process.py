import KratosMultiphysics
import KratosMultiphysics.SolidMechanicsApplication as KratosSolid

import assign_value_to_scalar_process as BaseProcess

def Factory(custom_settings, Model):
    if(type(custom_settings) != KratosMultiphysics.Parameters):
        raise Exception("Expected input shall be a Parameters object, encapsulating a json string")
    return AssignScalarToConditionsProcess(Model, custom_settings["Parameters"])

##all the processes python processes should be derived from "python_process"

class AssignScalarToConditionsProcess(BaseProcess.AssignValueToScalarProcess):
    def __init__(self, Model, custom_settings ):
        BaseProcess.AssignValueToScalarProcess.__init__(self, Model, custom_settings)        
        
    def ExecuteInitialize(self):

       if( self.interval_string == "initial" ):
            for cond in self.model_part.Conditions:
                cond.SetValue(self.var,self.value)

    def ExecuteInitializeSolutionStep(self):

        current_time = self.model_part.ProcessInfo[KratosMultiphysics.TIME]
        delta_time = self.model_part.ProcessInfo[KratosMultiphysics.DELTA_TIME]

        #arithmetic floating point tolerance
        tolerance = delta_time * 0.001;

        if( current_time >= (self.interval[0] - tolerance) and current_time <= (self.interval[1] + tolerance) ):
            if( self.function_string == "constant" ):
                for cond in self.model_part.Conditions:
                    cond.SetValue(self.var,self.value)
            else:        
                current_value = self.value * self.function(current_time)
                for cond in self.model_part.Conditions:
                    cond.SetValue(self.var,current_value)
    
        
    
