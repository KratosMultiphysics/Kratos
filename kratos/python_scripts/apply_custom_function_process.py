from KratosMultiphysics import *
import python_process
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
    def __init__(self, compiled_function, variable ):
        self.compiled_function = compiled_function
        self.variable = variable

    def f(self,x,y,z,t):
        return eval(self.compiled_function)
        
##all the processes python processes should be derived from "python_process"
class ApplyCustomFunctionProcess(python_process.PythonProcess):
    def __init__(self, model_part, mesh_id, variable_name, interval, is_fixed, function_string, free_outside_of_interval=True):
        python_process.PythonProcess.__init__(self) 
        
        self.model_part = model_part
        self.variable = globals().get(variable_name)
        self.mesh = model_part.GetMesh(mesh_id)
        #self.table = model_part.GetTable(table_id)
        self.interval = interval #expected to be a vector of size 2
        self.is_fixed = is_fixed
        self.free_outside_of_interval = free_outside_of_interval
        self.function_string = function_string
        
        ##compile function function_string
        self.compiled_function = compile(self.function_string, '', 'eval') #here we precompile the expression so that then it is much faster to execute it
        
        #construct a variable_utils object to speedup fixing
        self.variable_utils = VariableUtils()
        
        self.tmp = aux_object_cpp_callback(self.compiled_function, self.variable)
        self.cpp_apply_function_utility = PythonGenericFunctionUtility(self.variable, self.mesh.Nodes, self.tmp ) 
            
        print("finished construction of ApplyCustomFunctionProcess Process")
    
    
    def function(self, x,y,z,t):
        return eval(self.compiled_function) #this one is MUCH FASTER!
        #return eval(self.function_string)
      
            

    
    def ExecuteInitializeSolutionStep(self):
        current_time = self.model_part.ProcessInfo[TIME]
        
        if(current_time > self.interval[0] and  current_time<self.interval[1]):
            self.variable_utils.ApplyFixity(self.variable, self.is_fixed, self.mesh.Nodes)
            
            #create and pass around a list
            #tmp = aux_object(current_time, self.compiled_function, self.variable)   
            #self.variable_utils.ApplyVector(self.variable, list(map(tmp.freturn, self.mesh.Nodes)),  self.mesh.Nodes)
                
            #other option - also faster than the first for larger examples
            #tmp = aux_object(current_time, self.compiled_function, self.variable)   
            #any(map(tmp.f, self.mesh.Nodes))
        
            #faster for larger examples
            #for node in self.mesh.Nodes:
                #node.SetSolutionStepValue(self.variable,0,    self.function(node.X,node.Y,node.Z, current_time)    )
                
            #ultimate version. calling custom defined function from c++
            self.cpp_apply_function_utility.ApplyFunction(current_time)
            
    def ExecuteFinalizeSolutionStep(self):
        if self.free_outside_of_interval:
            current_time = self.model_part.ProcessInfo[TIME]
        
            if(current_time > self.interval[0] and  current_time<self.interval[1]):
                
                #here we free all of the nodes in the mesh
                fix = False
                self.variable_utils.ApplyFixity(self.variable, fix, self.mesh.Nodes)

    
    
    
def Factory(settings, Model):
    params = settings["parameters"]
    
    if(settings["process_name"] == "ApplyCustomFunctionProcess"):
        model_part = Model.get(  params.get( "model_part_name", "not found!!" ) , "model part not found" )
        mesh_id = params["mesh_id"]
        variable_name = params["variable_name"] 
        is_fixed = params["is_fixed"]
        #table_id = params["table_id"]
        interval = params["interval"]
        free_outside_of_interval = params["free_outside_of_interval"]
        function_string = params["f(x,y,z,t)="]
        
        if(len(interval) != 2):
            raise Exception( "interval size is expected to be two, while we got interval = "+str(interval));
        
        #here we could eventually sanitize the function
        
        return ApplyCustomFunctionProcess(model_part, mesh_id, variable_name, interval, is_fixed, function_string, free_outside_of_interval)
    else:
        raise Exception("trying to construct a Process with the wrong name!")