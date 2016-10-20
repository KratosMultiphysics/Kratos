import KratosMultiphysics
import KratosMultiphysics.SolidMechanicsApplication as KSM
import apply_vector_on_conditions_process

import math

def Factory(settings, Model):
    if(type(settings) != KratosMultiphysics.Parameters):
        raise Exception("Expected input shall be a Parameters object, encapsulating a json string")
    return ApplyVaryingVectorOverTimeIntervalProcess(Model, settings["Parameters"])

from apply_vector_on_conditions_process import ApplyVectorOnConditionsProcess as parent_process
from KratosMultiphysics import Process as base_process

##all the processes python processes should be derived from "python_process"
class ApplyVaryingVectorOverTimeIntervalProcess(parent_process):
    def __init__(self, Model, settings ):
        base_process.__init__(self)
##### 
#        -- SAMPLE --
#        "Parameters"            : {
#            "mesh_id"         : 0,
#            "model_part_name" : "LineLoad2D_bc_pressure",
#            "variable_name"   : "LINE_LOAD",
#            "direction"       : [0.0, -1.0, 0.0]
#            "factor1"         : 0.0,
#            "factor2"         : 0.2,
#            "interval"        : [0.0, 1.0]
#            "step_type"       : "linear"
#           }
#####

        self.Model = Model
 
        self.computation_model_part    = settings["model_part_name"]
        self.mesh_id        = settings["mesh_id"]
        self.variable_name  = settings["variable_name"]
                            
        self.curr_factor    = KratosMultiphysics.Parameters('{ "factor": 0.0 }');
        self.factor1        = settings["factor1"];
        self.factor2        = settings["factor2"];
        self.curr_direction = settings["direction"];
        self.interval       = settings["interval"]
        self.step_type      = settings["step_type"]
         
    def ExecuteInitialize(self):
        pass
 
    def ExecuteInitializeSolutionStep(self):
        current_time = self.Model[self.computation_model_part.GetString()].ProcessInfo[KratosMultiphysics.TIME]
        
        if( current_time < self.interval[0].GetDouble() ):   # before interval
            self.curr_factor = self.factor1
        elif( current_time > self.interval[1].GetDouble() ):   # after interval
            self.curr_factor = self.factor2
        else: # within interval
            if( self.step_type.GetString() == "constant" ):
                self.curr_factor = self.factor1
            
            elif( self.step_type.GetString() == "linear" ):
                self.curr_factor.SetDouble(
                                           self.factor1.GetDouble()
                                           + current_time * (
                                                             (self.factor2.GetDouble() - self.factor1.GetDouble())
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
        curr_step_params.AddValue("model_part_name",self.computation_model_part)
        curr_step_params.AddValue("mesh_id",self.mesh_id)
        curr_step_params.AddValue("variable_name",self.variable_name)
        curr_step_params.AddValue("modulus",self.curr_factor)
        curr_step_params.AddValue("direction",self.curr_direction)

        if( self.curr_factor.GetDouble() > 1e-9 ):
            parent_process.__init__(self, self.Model, curr_step_params)
            parent_process.ExecuteInitialize(self)
