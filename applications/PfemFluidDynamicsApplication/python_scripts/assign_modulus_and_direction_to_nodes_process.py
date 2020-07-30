from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7
#import kratos core and applications
import KratosMultiphysics

## This proces sets the value of a vector variable by a direction and a modulus.
## In this case, the fixicity is given set by deffault to true.
import math
import sys

from KratosMultiphysics.PfemFluidDynamicsApplication.assign_vector_components_to_nodes_process import AssignVectorComponentsToNodesProcess

def Factory(custom_settings, Model):
    if( not isinstance(custom_settings,KratosMultiphysics.Parameters) ):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return AssignModulusAndDirectionToNodesProcess(Model, custom_settings["Parameters"])

## All the processes python should be derived from "Process"
class AssignModulusAndDirectionToNodesProcess(AssignVectorComponentsToNodesProcess):
    def __init__(self, Model, custom_settings ):
        KratosMultiphysics.Process.__init__(self)

        ##settings string in json format
        default_settings = KratosMultiphysics.Parameters("""
        {
             "help" : "This process assigns a modulus and direction value to a vector variable",
             "model_part_name": "MODEL_PART_NAME",
             "variable_name": "VARIABLE_NAME",
             "modulus" : 0.0,
             "direction": [0.0, 0.0, 0.0],
             "compound_assignment": "direct",
             "constrained": false,
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
        if( not isinstance(self.var, KratosMultiphysics.Array1DVariable3) ):
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

        AssignVectorComponentsToNodesProcess.__init__(self, Model, params)
