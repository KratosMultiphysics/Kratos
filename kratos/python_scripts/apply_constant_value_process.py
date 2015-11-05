from KratosMultiphysics import *
import python_process

##all the processes python processes should be derived from "python_process"
class ApplyConstantScalarValue(python_process.PythonProcess):
    def __init__(self, model_part, variable_name, value, is_fixed, mesh_id=0 ):
        python_process.PythonProcess.__init__(self) 
        
        variable = globals().get(variable_name)

        for node in model_part.GetMesh(mesh_id).Nodes:
            if is_fixed:
                node.Fix(variable)
            
            node.SetSolutionStepValue(variable,0,value)
            
            
        print("finished construction of ApplyConstantScalarValue Process")
        
    def ExecuteInitialize(self):
        pass
    
    def ExecuteBeforeSolutionLoop(self):
        pass
    
    def ExecuteInitializeSolutionStep(self):
        pass
    
    def ExecuteFinalizeSolutionStep(self):
        pass
    
    def ExecuteBeforeOutputStep(self):
        pass
    
    def ExecuteAfterOutputStep(self):
        pass
        
    def ExecuteFinalize(self):
        pass
    
    def Clear(self):
        pass



class ApplyConstantVectorValue(python_process.PythonProcess):
    def __init__(self, model_part, variable_name, value, is_fixed_x, is_fixed_y, is_fixed_z, mesh_id=0):
        python_process.PythonProcess.__init__(self) 

        #print("value = ", value, len(value))
        variable = globals().get(variable_name)
        varx = globals().get(variable_name+"_X")
        vary = globals().get(variable_name+"_Y")
        varz = globals().get(variable_name+"_Z")
        
        if(len(value) != 3):
            raise Exception("sorry the value to be applied should be a vector of size 3. Instead we got", value)
        
        for node in model_part.GetMesh(mesh_id).Nodes:
            if is_fixed_x:
                node.Fix(varx)
            if is_fixed_y:
                node.Fix(vary)
            if is_fixed_z:
                node.Fix(varz)
                
            node.SetSolutionStepValue(variable,0,value)
        print("finished construction of ApplyConstantVectorValue Process")
        
    def ExecuteInitialize(self):
        pass
    
    def ExecuteBeforeSolutionLoop(self):
        pass
    
    def ExecuteInitializeSolutionStep(self):
        pass
    
    def ExecuteFinalizeSolutionStep(self):
        pass
    
    def ExecuteBeforeOutputStep(self):
        pass
    
    def ExecuteAfterOutputStep(self):
        pass
        
    def ExecuteFinalize(self):
        pass
    
    def Clear(self):
        pass
    
    
def Factory(settings, Model):
    params = settings["parameters"]
    model_part = Model.get(  params.get( "model_part_name", "not found!!" ) , "model part not found" )
    mesh_id = settings["mesh_id"]
    
    if(settings["process_name"] == "ApplyConstantScalarValue"):
        variable_name = params["variable_name"] 
        value = params["value"]
        is_fixed = params["is_fixed"]
        
        return ApplyConstantScalarValue(model_part, variable_name, value,is_fixed)
    elif(settings["process_name"] == "ApplyConstantVectorValue"):
        variable_name = params["variable_name"]  
        value = params["value"]
        #print("value is ", value,len(value))
        
        is_fixed_x = params["is_fixed_x"]
        is_fixed_y = params["is_fixed_y"]
        is_fixed_z = params["is_fixed_z"]
        
        return ApplyConstantVectorValue(model_part, variable_name, value,is_fixed_x, is_fixed_y, is_fixed_z)