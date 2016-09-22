import KratosMultiphysics
import impose_vector_value_by_components_process

## This proces sets the value of a vector variable component-by-component.
## In this case, the fixicity is given by the user and some of the components may not be fixed.

def Factory(settings, Model):
    if(type(settings) != KratosMultiphysics.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return ImposeVectorValueByComponentsOverTimeProcess(Model, settings["Parameters"])

from impose_vector_value_by_components_process import ImposeVectorValueByComponentsProcess as parent_process
from KratosMultiphysics import Process as base_process

## All the processes python processes should be derived from "python_process"
class ImposeVectorValueByComponentsOverTimeProcess(parent_process):
    def __init__(self, Model, settings ):
        base_process.__init__(self)
##### 
#        -- SAMPLE --
#       "Parameters"            : {
#            "model_part_name" : "el_model_part",
#            "mesh_id"         : 0,
#            "variable_name"   : "DISPLACEMENT",
#            "is_fixed_x"      : false,
#            "is_fixed_y"      : true,
#            "is_fixed_z"      : false,
#            "value1"          : [0.0,-0.14,0.0],
#            "value2"          : [0.0,-0.14,0.0],
#            "interval"        : [0.0, 1.0],
#            "step_type"       : "linear"
#           }
#####

        self.Model = Model

        self.model_part    = settings["model_part_name"]
        self.mesh_id       = settings["mesh_id"]
        self.variable_name = settings["variable_name"]
        self.is_fixed_x    = settings["is_fixed_x"]
        self.is_fixed_y    = settings["is_fixed_y"]
        self.is_fixed_z    = settings["is_fixed_z"]
    
        self.value1        = settings["value1"]
        print( self.value1.PrettyPrintJsonString() )
        self.value2        = settings["value2"]
        print( self.value2.PrettyPrintJsonString() )
        self.interval      = settings["interval"]
        self.step_type     = settings["step_type"]
       
    def ExecuteInitialize(self):
        pass
        
    def ExecuteInitializeSolutionStep(self):
        model_part = self.Model[self.model_part.GetString()]
        current_time = model_part.ProcessInfo[KratosMultiphysics.TIME]
        
        curr_value = KratosMultiphysics.Parameters('{ "value": [0.0, 0.0, 0.0] }')
        
        if( current_time < self.interval[0].GetDouble() ):   # before interval
            print( "\033[93m :: Taking initial value :: \033[0m" )
            for i in range(2):
                curr_value["value"][i].SetDouble( self.value1[0].GetDouble() )
                
        elif( current_time > self.interval[1].GetDouble() ): # after interval
            print( "\033[93m :: Taking final value :: \033[0m" )
            for i in range(2):
                curr_value["value"][i].SetDouble( self.value2[0].GetDouble() )

        else:                                                # within interval
            if( self.step_type.GetString() == "constant" ):
                print( "\033[93m :: CONSTANT: Taking initial value :: \033[0m" )
                for i in range(2):
                    curr_value["value"][i].SetDouble( self.value1[0].GetDouble() )
            
            elif( self.step_type.GetString() == "linear" ):
                print( "\033[93m :: LINEAR: Interpolating current value :: \033[0m" )
                for i in range(2):
                    curr_value["value"][i].SetDouble(
                                                     self.value1[i].GetDouble()
                                                     + current_time * (
                                                                       (self.value2[i].GetDouble() - self.value1[i].GetDouble())
                                                                       / 
                                                                       (self.interval[1].GetDouble() - self.interval[0].GetDouble())
                                                                       )
                                                     )
            else:
                raise Exception("Only constant and linear variation is implemented")
        
        print( "#####################################################" )
        print( "\033[92m", "::: Time: ", current_time, " :::", "\033[0m" )
        print( "#####################################################" )
        
        curr_step_params = KratosMultiphysics.Parameters("{}")
        curr_step_params.AddValue("model_part_name",self.model_part)
        curr_step_params.AddValue("mesh_id",self.mesh_id)
        curr_step_params.AddValue("variable_name",self.variable_name)
        curr_step_params.AddValue("is_fixed_x",self.is_fixed_x)
        curr_step_params.AddValue("is_fixed_y",self.is_fixed_y)
        curr_step_params.AddValue("is_fixed_z",self.is_fixed_z)
        curr_step_params.AddValue("value",curr_value["value"])

        parent_process.__init__(self, self.Model, curr_step_params)
        parent_process.ExecuteInitialize(self)
