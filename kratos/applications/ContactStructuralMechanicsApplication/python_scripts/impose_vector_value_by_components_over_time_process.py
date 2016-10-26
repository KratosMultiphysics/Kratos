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
#            "model_part_name"  : "el_model_part",
#            "mesh_id"          : 0,
#            "variable_name"    : "DISPLACEMENT",
#            "vectors"          : [ [0.0,-0.14,0.0],[0.0,-0.14,0.0] ],
#            "intervals"        : [0.0, 1.0],
#            "step_type"        : "linear"
#           }
#####

        self.Model = Model

        self.model_part    = settings["model_part_name"]
        self.mesh_id       = settings["mesh_id"]
        self.variable_name = settings["variable_name"]
    
        self.vectors   = settings["vectors"]
        self.intervals = settings["intervals"]
        self.step_type = settings["step_type"]
        
        # checks
        assert( self.vectors.size( ) == self.intervals.size( ) )
        assert( self.intervals[self.intervals.size( )-1].GetDouble() >= 0.0 )
       
    def ExecuteInitialize(self):
        pass
        
    def ExecuteInitializeSolutionStep(self):
        current_time = self.Model[self.model_part.GetString()].ProcessInfo[KratosMultiphysics.TIME]
        
        curr_value = KratosMultiphysics.Parameters('{ "value": [0.0, 0.0, 0.0] }')
        curr_is_fixed_x = KratosMultiphysics.Parameters('{ "is_fixed_x": false }')
        curr_is_fixed_y = KratosMultiphysics.Parameters('{ "is_fixed_y": false }')
        curr_is_fixed_z = KratosMultiphysics.Parameters('{ "is_fixed_z": false }')
        
        n_intervals = self.intervals.size( ); # number of steps
        
        # outside interval
        if( current_time < self.intervals[0].GetDouble() or
            current_time > self.intervals[n_intervals-1].GetDouble() ):
            curr_is_fixed_x.SetBool(False)
            curr_is_fixed_y.SetBool(False)
            curr_is_fixed_z.SetBool(False)
                
        # inside interval
        else:
            # fix the DOFs
            curr_is_fixed_x.SetBool(True)
            curr_is_fixed_y.SetBool(True)
            curr_is_fixed_z.SetBool(True)
            
            # interpolate the value
            for n in range(1, n_intervals):
                if (self.intervals[n].GetDouble() > current_time):
                    break;
            
            for i in range(2):
                if( abs( self.vectors[n][i].GetDouble() - self.vectors[n-1][i].GetDouble() ) < 1e-6 ):
                    curr_value["value"][i].SetDouble( self.vectors[n-1][i].GetDouble() )
                else:
                    y1 = self.vectors[n-1][i].GetDouble();
                    y2 = self.vectors[n][i].GetDouble();
                    dy = y2 - y1;
                    t1 = self.intervals[n-1].GetDouble();
                    t2 = self.intervals[n].GetDouble();
                    dt = t2 - t1;
                    t = current_time;
                    
                    if( self.step_type.GetString() == "linear" ):
                        curr_value["value"][i].SetDouble( y1  + ( t - t1 ) * ( dy/dt ) )
                        
                    elif( self.step_type.GetString() == "smooth" ):
                        from numpy import tanh as tanh  # we use hyperbolic tan to define a smooth step
                        t_av = 0.5 * (t1 + t2);
                        y_av = 0.5 * (y1 + y2);
                        c_slope = 0.15;
                        curr_value["value"][i].SetDouble( y_av + 0.5 * dy * tanh( (t-t_av)/(c_slope*dt) ) )
                        
                    else:
                        raise Exception("Only linear and smooth are valid inputs")
                
        print( "#####################################################" )
        print( "\033[92m" )
        print( "::: Time: ", current_time, " :::" )
        print( curr_value["value"].PrettyPrintJsonString() )
        print( "t_{n-1} = ", self.intervals[n-1].PrettyPrintJsonString() )
        print( "t_{n}   = ", self.intervals[n  ].PrettyPrintJsonString() )
        print( "v_{n-1} = ", self.vectors[n-1].PrettyPrintJsonString() )
        print( "v_{n}   = ", self.vectors[n  ].PrettyPrintJsonString() )
        print( "\033[0m" )
        print( "#####################################################" )
        
        curr_step_params = KratosMultiphysics.Parameters("{}")
        curr_step_params.AddValue("model_part_name",self.model_part)
        curr_step_params.AddValue("mesh_id",self.mesh_id)
        curr_step_params.AddValue("variable_name",self.variable_name)
        curr_step_params.AddValue("is_fixed_x",curr_is_fixed_x)
        curr_step_params.AddValue("is_fixed_y",curr_is_fixed_y)
        curr_step_params.AddValue("is_fixed_z",curr_is_fixed_z)
        curr_step_params.AddValue("value",curr_value["value"])

        parent_process.__init__(self, self.Model, curr_step_params)
        parent_process.ExecuteInitialize(self)
