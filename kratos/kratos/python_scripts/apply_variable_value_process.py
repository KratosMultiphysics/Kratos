from KratosMultiphysics import *
import python_process

##all the processes python processes should be derived from "python_process"
class ApplyVariableScalarValue(python_process.PythonProcess):
    def __init__(self, model_part, variable_name, table_id, interval, is_fixed, free_outside_of_interval=True, mesh_id=0 ):
        python_process.PythonProcess.__init__(self) 
        
        self.model_part = model_part
        self.variable = globals().get(variable_name)
        self.mesh = model_part.GetMesh(mesh_id)
        self.table = model_part.GetTable(table_id)
        self.interval = interval #expected to be a vector of size 2
        self.is_fixed = is_fixed
        self.free_outside_of_interval = free_outside_of_interval
            
        print("finished construction of ApplyVariableScalarValue Process")
            
    def ExecuteInitializeSolutionStep(self):
        current_time = self.model_part.ProcessInfo[TIME]
        
        if(current_time > self.interval[0] and  current_time<self.interval[1]):
            if self.is_fixed:
                for node in self.mesh.Nodes:
                    node.Fix(self.variable)
            
            value = self.table.GetValue(current_time)
            
            for node in self.mesh.Nodes:
                node.SetSolutionStepValue(self.variable,0,value)
            
    def ExecuteFinalizeSolutionStep(self):
        if self.free_outside_of_interval:
            current_time = self.model_part.ProcessInfo[TIME]
        
            if(current_time > self.interval[0] and  current_time<self.interval[1]):
                for node in self.mesh.Nodes:
                    node.Free(self.variable)
    
    
    
def Factory(settings, Model):
    params = settings["parameters"]
    
    if(settings["process_name"] == "ApplyVariableScalarValue"):
        model_part = Model.get(  params.get( "model_part_name", "not found!!" ) , "model part not found" )
        mesh_id = settings["mesh_id"]
        variable_name = params["variable_name"] 
        value = params["value"]
        is_fixed = params["is_fixed"]
        table_id = params["table_id"]
        interval = params["interval"]
        free_outside_of_interval = params["free_outside_of_interval"]
        
        return ApplyVariableScalarValue(model_part, variable_name, table_id, interval, is_fixed, free_outside_of_interval, mesh_id )
    else:
        raise Exception("trying to construct a Process with the wrong name!")