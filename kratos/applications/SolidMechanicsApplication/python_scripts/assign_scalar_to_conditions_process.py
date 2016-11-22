import KratosMultiphysics
import KratosMultiphysics.SolidMechanicsApplication as KratosSolid

## This proces sets the value of a scalar variable to conditions

import assign_scalar_to_nodes_process as BaseProcess

def Factory(custom_settings, Model):
    if(type(custom_settings) != KratosMultiphysics.Parameters):
        raise Exception("Expected input shall be a Parameters object, encapsulating a json string")
    return AssignScalarToConditionsProcess(Model, custom_settings["Parameters"])


class AssignScalarToConditionsProcess(BaseProcess.AssignScalarToNodesProcess):
    def __init__(self, Model, custom_settings ):
        BaseProcess.AssignScalarToNodesProcess.__init__(self, Model, custom_settings)        
         
    def ExecuteInitialize(self):

        # set processes
        params = KratosMultiphysics.Parameters("{}")           
        params.AddValue("model_part_name", self.settings["model_part_name"])
        params.AddValue("mesh_id", self.settings["mesh_id"])
        params.AddValue("variable_name", self.settings["variable_name"])
        
        if( self.value_is_numeric ):
            print(" numeric value ",self.value)
            params.AddValue("value", self.settings["value"])
            self.AssignValueProcess = KratosSolid.AssignScalarToConditionsProcess(self.model_part, params)
        else:
            print(" function value ", self.value)
            self.AssignValueProcess = KratosSolid.AssignScalarFieldToConditionsProcess(self.model_part, self.compiled_function, "function", self.value_is_spatial_function, params)
            
        if( self.IsInsideInterval() and self.interval_string == "initial" ):
            self.AssignValueProcess.Execute()
                    
    def ExecuteInitializeSolutionStep(self):

        if self.IsInsideInterval():
            self.AssignValueProcess.Execute()

    def ExecuteFinalizeSolutionStep(self):
        pass
            

    
