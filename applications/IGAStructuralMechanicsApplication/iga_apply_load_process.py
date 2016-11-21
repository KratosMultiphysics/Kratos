import KratosMultiphysics
import python_process

##all the processes python processes should be derived from "python_process"
class IGAApplyLoad(python_process.PythonProcess):
    def __init__(self, model_part, variable_name, factor, direction, mesh_id=0 ):
        python_process.PythonProcess.__init__(self) 
        
        variable = globals().get(variable_name)
        #print(model_part)
        for condition in model_part.GetMesh(mesh_id).Conditions:
            condition.SetValue(KratosMultiphysics.SolidMechanicsApplication.POINT_LOAD,[direction[0].GetDouble(),direction[1].GetDouble(),direction[2].GetDouble()])
        #print(direction[0].GetDouble(),direction[1].GetDouble(),direction[2].GetDouble())
        print("Finished construction of IGAApplyLoad Process")
        
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
	print(params["model_part_name"].GetString())
	model_part = Model.get(  params["model_part_name"].GetString() , "model part not found" )
	mesh_id = params["mesh_id"]

	#if(settings["process_name"] == "IGAApplyLoad"):
	variable_name = params["variable_name"] 
	factor = params["factor"]
	direction = params["direction"]

	return IGAApplyLoad(model_part, variable_name, factor, direction)