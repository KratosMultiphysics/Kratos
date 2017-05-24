import KratosMultiphysics
import impose_scalar_value_process

## This proces sets the value of a vector variable component-by-component.
## In this case, the fixicity is given by the user and some of the components may not be fixed.

def Factory(settings, Model):
    if(type(settings) != KratosMultiphysics.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return ImposeScalarValueOverTimeProcess(Model, settings["Parameters"])

from impose_scalar_value_process import ImposeScalarValueProcess as parent_process
from KratosMultiphysics import Process as base_process

## All the processes python processes should be derived from "python_process"
class ImposeScalarValueOverTimeProcess(parent_process):
    def __init__(self, Model, settings ):
        base_process.__init__(self)
##### 
#        -- SAMPLE --
#       "Parameters"            : {
#            "model_part_name"  : "el_model_part",
#            "mesh_id"          : 0,
#            "variable_name"    : "DISPLACEMENT_X",
#            "values"          :  [0.0,-0.14,0.0 ],
#            "intervals"        : [0.0, 1.0],
#            "step_type"        : "linear"
#           }
#####

        self.Model = Model

        self.model_part    = settings["model_part_name"]
        self.mesh_id       = settings["mesh_id"]
        self.variable_name = settings["variable_name"]
    
        self.values   = settings["values"]
        self.intervals = settings["intervals"]
        self.step_type = settings["step_type"]
        
        # checks
        assert( self.values.size( ) == self.intervals.size( ) )
        assert( self.intervals[self.intervals.size( )-1].GetDouble() >= 0.0 )
       
    def ExecuteInitialize(self):
        pass
        
    def ExecuteInitializeSolutionStep(self):
        current_time = self.Model[self.model_part.GetString()].ProcessInfo[KratosMultiphysics.TIME]
        
        curr_value = KratosMultiphysics.Parameters('{ "value": 0.0 }')
        curr_is_fixed = KratosMultiphysics.Parameters('{ "is_fixed": false }')
        n_intervals = self.intervals.size( ); # number of steps
        
        # outside interval
        if( current_time < self.intervals[0].GetDouble() or
            current_time > self.intervals[n_intervals-1].GetDouble() ):
            curr_is_fixed.SetBool(False)
                
        # inside interval
        else:
            # fix the DOFs
            curr_is_fixed.SetBool(True)
            
            # interpolate the value
            for n in range(1, n_intervals):
                if (self.intervals[n].GetDouble() > current_time):
                    break;
            
            if( abs( self.values[n].GetDouble() - self.values[n-1].GetDouble() ) < 1e-6 ):
                curr_value["value"].SetDouble( self.values[n-1].GetDouble() )
            else:
                y1 = self.values[n-1].GetDouble();
                y2 = self.values[n].GetDouble();
                dy = y2 - y1;
                t1 = self.intervals[n-1].GetDouble();
                t2 = self.intervals[n].GetDouble();
                dt = t2 - t1;
                t = current_time;
                
                if( self.step_type.GetString() == "linear" ):
                    curr_value["value"].SetDouble( y1  + ( t - t1 ) * ( dy/dt ) )
                    
                elif( self.step_type.GetString() == "smooth" ):
                    from math import tanh as tanh  # we use hyperbolic tan to define a smooth step
                    t_av = 0.5 * (t1 + t2);
                    y_av = 0.5 * (y1 + y2);
                    c_slope = 0.15;
                    curr_value["value"].SetDouble( y_av + 0.5 * dy * tanh( (t-t_av)/(c_slope*dt) ) )
                    
                else:
                    raise Exception("Only linear and smooth are valid inputs")
                
        print( "#####################################################" )
        print( "\033[92m" )
        print( "::: Time: ", current_time, " :::" )
        print( curr_value["value"].PrettyPrintJsonString() )
        print( "t_{n-1} = ", self.intervals[n-1].PrettyPrintJsonString() )
        print( "t_{n}   = ", self.intervals[n  ].PrettyPrintJsonString() )
        print( "v_{n-1} = ", self.values[n-1].PrettyPrintJsonString() )
        print( "v_{n}   = ", self.values[n  ].PrettyPrintJsonString() )
        print( "\033[0m" )
        print( "#####################################################" )
        
        curr_step_params = KratosMultiphysics.Parameters("{}")
        curr_step_params.AddValue("model_part_name",self.model_part)
        curr_step_params.AddValue("mesh_id",self.mesh_id)
        curr_step_params.AddValue("variable_name",self.variable_name)
        curr_step_params.AddValue("is_fixed",curr_is_fixed)
        curr_step_params.AddValue("value",curr_value["value"])

        parent_process.__init__(self, self.Model, curr_step_params)
        parent_process.ExecuteInitialize(self)
