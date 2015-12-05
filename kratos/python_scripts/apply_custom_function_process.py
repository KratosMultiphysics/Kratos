from KratosMultiphysics import *
import python_process
from math import *

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
            
        print("finished construction of ApplyCustomFunctionProcess Process")
        
    def function(self, x,y,z,t):
        import math
        value = eval(self.function_string)
        return value
            
    def ExecuteInitializeSolutionStep(self):
        current_time = self.model_part.ProcessInfo[TIME]
        
        if(current_time > self.interval[0] and  current_time<self.interval[1]):
            if self.is_fixed:
                for node in self.mesh.Nodes:
                    node.Fix(self.variable)
                       
            for node in self.mesh.Nodes:
                value = self.function(node.X,node.Y,node.Z, current_time)
                
                node.SetSolutionStepValue(self.variable,0,value)
            
    def ExecuteFinalizeSolutionStep(self):
        if self.free_outside_of_interval:
            current_time = self.model_part.ProcessInfo[TIME]
        
            if(current_time > self.interval[0] and  current_time<self.interval[1]):
                for node in self.mesh.Nodes:
                    node.Free(self.variable)
    
    
    
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