import KratosMultiphysics 
import sys
from math import *

class aux_object:
    def __init__(self, t, compiled_function, variable ):
        self.t = t
        self.compiled_function = compiled_function
        self.variable = variable
    
    def f(self,node):
        x = node.X
        y = node.Y
        z = node.Z
        t = self.t
        node.SetSolutionStepValue(self.variable,0, eval(self.compiled_function) )

    def freturn(self,node):
        x = node.X
        y = node.Y
        z = node.Z
        t = self.t
        return eval(self.compiled_function)

class aux_object_cpp_callback:
    def __init__(self, compiled_function ):
        self.compiled_function = compiled_function

    def f(self,x,y,z,t):
        return eval(self.compiled_function)
        #return 0
        
def Factory(settings, Model):
    if(type(settings) != KratosMultiphysics.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return ApplyCustomFunctionProcess(Model, settings["Parameters"])

    
##all the processes python processes should be derived from "python_process"
class ApplyCustomFunctionProcess(KratosMultiphysics.Process):
    def __init__(self, Model, settings ):

        default_settings = KratosMultiphysics.Parameters("""
            {
                "mesh_id"         : 0,
                "model_part_name" : "Outlet2D_Outlet_pressure_Auto1",
                "variable_name"   : "TEMPERATURE",
                "interval"        : [0.0, 10.0],
                "is_fixed"		  : true,
                "free_outside_of_interval" : false,
                "f(x,y,z,t)"      : "x+y"
            }
            """
            )
            
        settings.ValidateAndAssignDefaults(default_settings)

        KratosMultiphysics.Process.__init__(self) 
        self.model_part = Model[settings["model_part_name"].GetString()]
        self.variable = getattr(KratosMultiphysics, settings["variable_name"].GetString())
        self.mesh = self.model_part.GetMesh(settings["mesh_id"].GetInt())
        #self.table = model_part.GetTable(table_id)
        self.interval = settings["interval"] #expected to be a vector of size 2
        self.is_fixed = settings["is_fixed"].GetBool()
        self.free_outside_of_interval = settings["free_outside_of_interval"].GetBool()
        self.function_string = settings["f(x,y,z,t)"].GetString()
        
        ##compile function function_string
        if (sys.version_info > (3, 0)):
            self.compiled_function = compile(self.function_string, '', 'eval', optimize=2) #here we precompile the expression so that then it is much faster to execute it
        else:
            self.compiled_function = compile(self.function_string, '', 'eval')
        
        #construct a variable_utils object to speedup fixing
        self.variable_utils = KratosMultiphysics.VariableUtils()
        
        self.tmp = KratosMultiphysics.PythonGenericFunctionUtility(aux_object_cpp_callback(self.compiled_function))
        self.cpp_apply_function_utility = KratosMultiphysics.ApplyFunctionToNodesUtility(self.mesh.Nodes, self.tmp ) 
            
        print("finished construction of ApplyCustomFunctionProcess Process")
    
    def function(self, x,y,z,t):
        return eval(self.compiled_function) #this one is MUCH FASTER!
        #return eval(self.function_string)
    
    def ExecuteInitializeSolutionStep(self):
        current_time = self.model_part.ProcessInfo[KratosMultiphysics.TIME]
        
        if(current_time > self.interval[0].GetDouble() and  current_time<self.interval[1].GetDouble()):
            self.variable_utils.ApplyFixity(self.variable, self.is_fixed, self.mesh.Nodes)
            
            #create and pass around a list
            #tmp = aux_object(current_time, self.compiled_function, self.variable)   
            #self.variable_utils.ApplyVector(self.variable, Vector( list(map(tmp.freturn, self.mesh.Nodes)) ),  self.mesh.Nodes)
                
            #other option - also faster than the first for larger examples
            #tmp = aux_object(current_time, self.compiled_function, self.variable)   
            #any(map(tmp.f, self.mesh.Nodes))
        
            #faster for larger examples
            #for node in self.mesh.Nodes:
                #node.SetSolutionStepValue(self.variable,0,    self.function(node.X,node.Y,node.Z, current_time)    )
                
            #ultimate version. calling custom defined function from c++
            self.cpp_apply_function_utility.ApplyFunction(self.variable, current_time)
            
    def ExecuteFinalizeSolutionStep(self):
        if self.free_outside_of_interval:
            current_time = self.model_part.ProcessInfo[KratosMultiphysics.TIME]
        
            if(current_time > self.interval[0].GetDouble() and  current_time<self.interval[1].GetDouble()):
                
                #here we free all of the nodes in the mesh
                fix = False
                self.variable_utils.ApplyFixity(self.variable, fix, self.mesh.Nodes)
    
# def Factory(settings, Model):
    # params = settings["parameters"]
    
    # if(settings["process_name"].GetString() == "ApplyCustomFunctionProcess"):
        # model_part = Model[params["model_part_name"].GetString()]
        # mesh_id = params["mesh_id"].GetInt()
        # variable_name = params["variable_name"].GetString() 
        # is_fixed = params["is_fixed"].GetBool()
        #table_id = params["table_id"]
        # interval = [params["interval"][0].GetDouble(), params["interval"][1].GetDouble()]
        # free_outside_of_interval = params["free_outside_of_interval"].GetBool()
        # function_string = params["f(x,y,z,t)="].GetString() 
        
        # if(len(interval) != 2):
            # raise Exception( "interval size is expected to be two, while we got interval = "+str(interval));
        
        #here we could eventually sanitize the function
        
        # return ApplyCustomFunctionProcess(model_part, mesh_id, variable_name, interval, is_fixed, function_string, free_outside_of_interval)
    # else:
        # raise Exception("trying to construct a Process with the wrong name!")