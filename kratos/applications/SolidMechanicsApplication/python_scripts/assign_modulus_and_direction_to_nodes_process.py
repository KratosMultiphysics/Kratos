import KratosMultiphysics
import KratosMultiphysics.SolidMechanicsApplication as KratosSolid

## This proces sets the value of a vector variable by a direction and a modulus.
## In this case, the fixicity is given set by deffault to true.
import math
import sys

import assign_vector_components_to_nodes_process as BaseProcess

def Factory(custom_settings, Model):
    if(type(custom_settings) != KratosMultiphysics.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return AssignModulusAndDirectionToNodesProcess(Model, custom_settings["Parameters"])

## All the processes python processes should be derived from "python_process"
class AssignModulusAndDirectionToNodesProcess(BaseProcess.AssignVectorComponentsToNodesProcess):
    def __init__(self, Model, custom_settings ):
        KratosMultiphysics.Process.__init__(self)
        
        ##settings string in json format
        default_settings = KratosMultiphysics.Parameters("""
        {
             "help" : "This process assigns a modulus and direction value to a vector variable",
             "model_part_name": "MODEL_PART_NAME",
             "mesh_id": 0,
             "variable_name": "VARIABLE_NAME",
             "modulus" : 0.0,
             "direction": [0.0, 0.0, 0.0],
             "constrained": true,
             "interval": [0.0, "End"],
             "local_axes" : {}
        }
        """)

        #trick to allow "value" to be a string or a double value
        if(custom_settings.Has("modulus")):
            if(custom_settings["modulus"].IsString()):
                default_settings["modulus"].SetString("0.0")
                
        ##overwrite the default settings with user-provided parameters
        self.settings = custom_settings
        self.settings.ValidateAndAssignDefaults(default_settings)
                
        ##check if variable type is a vector
        self.var = KratosMultiphysics.KratosGlobals.GetVariable(self.settings["variable_name"].GetString())
        if( type(self.var) != KratosMultiphysics.Array1DVariable3 ):
            raise Exception("Variable type is incorrect. Must be a three-component vector.")

        ##check normalized direction
        direction   = []
        scalar_prod = 0 
        for i in range(0, self.settings["direction"].size() ):
            direction.append( self.settings["direction"][i].GetDouble() )
            scalar_prod = scalar_prod + direction[i]*direction[i]
            
        norm = math.sqrt(scalar_prod)
        
        self.value = []
        if( norm != 0.0 ):
            for j in direction:
                self.value.append( j/norm )
        else:
            for j in direction:
                self.value.append(0.0)

        ## set the value
        self.value_is_numeric = False

        if self.settings["modulus"].IsNumber():
            self.value_is_numeric = True
            
            modulus = self.settings["modulus"].GetDouble()
            for i in range(0, len(self.value)):
                self.value[i] *= modulus

        else:            
            self.function_expression = self.settings["modulus"].GetString()
            
            counter = 0
            for i in self.value:
                self.value[counter] = str(i) + "*(" + self.function_expression +")" 
                counter+=1
 
                    
        ###set vector components process
        params = KratosMultiphysics.Parameters("{}")           
        params.AddValue("model_part_name", self.settings["model_part_name"])
        params.AddValue("mesh_id", self.settings["mesh_id"])
        params.AddValue("variable_name", self.settings["variable_name"])

        params.AddEmptyValue("value")
        params.__setitem__("value", self.settings["direction"])
        counter = 0
        for i in self.value:
            if( self.value_is_numeric ):
                params["value"][counter].SetDouble(i)
            else:
                params["value"][counter].SetString(i)
            counter+=1

        params.AddValue("constrained", self.settings["constrained"])
        params.AddValue("interval",self.settings["interval"])
        params.AddValue("local_axes", self.settings["local_axes"])

        BaseProcess.AssignVectorComponentsToNodesProcess.__init__(self, Model, params)
        

            

        
        
        
