from KratosMultiphysics import *
import python_process

##all the processes python processes should be derived from "python_process"
class ApplyAcceleration(python_process.PythonProcess):
    def __init__(self, model_part, acceleration_variable_name, table_id, interval, is_fixed, free_outside_of_interval=True, mesh_id=0 ):
        python_process.PythonProcess.__init__(self) 
      
        if(acceleration_variable_name == "ACCELERATION_X"):
            variable = DISPLACEMENT_X
        if(acceleration_variable_name == "ACCELERATION_Y"):
            variable = DISPLACEMENT_Y
        if(acceleration_variable_name == "ACCELERATION_Z"):
            variable = DISPLACEMENT_Z
        
        #here compute a displacement table by integrating the acceleration table
        #TODO: implement this!!
        #then save it to a new displacement_table_id
        
        #now construct an internal process to apply the integrated displacement
        import apply_variable_value_process
        self.aux_process = apply_variable_value_process.ApplyVariableScalarValue(model_part, variable, displacement_table_id, interval, is_fixed, free_outside_of_interval mesh_id )
      
        print("finished construction of ApplyAcceleration Process")
            
    def ExecuteInitializeSolutionStep(self):
        self.aux_process.ExecuteInitializeSolutionStep(self)

            
    def ExecuteFinalizeSolutionStep(self):
        self.aux_process.ExecuteFinalizeSolutionStep(self)
    
    
    
def Factory(settings, Model):
    params = settings["parameters"]
    
    if(settings["process_name"] == "ApplyAcceleration"):
        model_part = Model.get(  params.get( "model_part_name", "not found!!" ) , "model part not found" )
        mesh_id = settings["mesh_id"]
        acceleration_variable_name = params["acceleration_variable_name"] 
        value = params["value"]
        is_fixed = params["is_fixed"]
        table_id = params["table_id"]
        interval = params["interval"]
        free_outside_of_interval = params["free_outside_of_interval"]
        
        return ApplyVariableScalarValue(model_part, variable_name, table_id, interval, is_fixed, free_outside_of_interval, mesh_id )
    else:
        raise Exception("trying to construct a Process with the wrong name (should be ApplyAcceleration)!")