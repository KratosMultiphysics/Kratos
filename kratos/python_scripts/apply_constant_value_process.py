from KratosMultiphysics import *
import python_process

##all the processes python processes should be derived from "python_process"
class ApplyConstantScalarValue(python_process.PythonProcess):
    def __init__(self, variable, value, is_fixed):
        python_process.PythonProcess.__init__(self) 
        self.is_fixed = is_fixed
        self.value = value
        self.variable = variable
        
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
        print("inside apply_constant_scalar_value_process")
        
    def ExecuteFinalize(self):
        pass
    
    def Clear(self):
        pass



class ApplyConstantVectorValue(python_process.PythonProcess):
    def __init__(self, variable, value, is_fixed_x, is_fixed_y, is_fixed_z):
        python_process.PythonProcess.__init__(self) 
        self.is_fixed_x = is_fixed_x
        self.is_fixed_y = is_fixed_y
        self.is_fixed_z = is_fixed_z
        self.value = value
        self.variable = variable
        
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
        print("inside apply_constant_scalar_value_process")
        
    def ExecuteFinalize(self):
        pass
    
    def Clear(self):
        pass
    
    
def Factory(settings, Model):
    #get the mesh and model part of interest
    model_part = globals()[ Model.get( "model_part_name", "not found!!" ) ] 
    mesh_id = settings["mesh_id"]
    
    if(settings["process_name"] == "ApplyConstantScalarValue"):
        #query for the list of actual parameters to be employed in construction
        params = settings["parameters"]
        
        variable = globals()[ params["variable_name"] ] 
        value = params["value"]
        is_fixed = settings["is_fixed"]
        
        return ApplyConstantScalarValue(variable, value,is_fixed)
    elif(settings["process_name"] == "ApplyConstantVectorValue"):
        #query for the list of actual parameters to be employed in construction
        params = settings["parameters"]
        
        variable = globals()[ params["variable_name"] ] 
        value = params["value"]
        if(len(value) != 3):
            raise Exception("value should be a vector of size 3. Instead we got ", value)
        is_fixed_x = settings["is_fixed_x"]
        is_fixed_y = settings["is_fixed_y"]
        is_fixed_z = settings["is_fixed_z"]
        
        return ApplyConstantVectorValue(variable, value,is_fixed_x, is_fixed_y, is_fixed_z)