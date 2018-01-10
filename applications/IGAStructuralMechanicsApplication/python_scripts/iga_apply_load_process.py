import KratosMultiphysics
import python_process

##all the processes python processes should be derived from "python_process"
class IGAApplyLoad(python_process.PythonProcess):
    def __init__(self, model_part, variables, mesh_id=0 ):
        python_process.PythonProcess.__init__(self) 
        
        for condition in model_part.Conditions:
            for variable_key in variables:
                condition.SetValue(eval(variable_key), variables[variable_key])


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
	#print(params["model_part_name"].GetString())
	#model_part = Model.get(  params["model_part_name"].GetString() , "model part not found" )
	model_part = Model[params["model_part_name"].GetString()].GetSubModelPart(params["sub_model_part_name"].GetString())# , "model part not found" )
	mesh_id = params["mesh_id"]

	variables = {}
	#if(settings["process_name"] == "IGAApplyLoad"):
	for variable_i in range (0,params["variables"].size()):
		variable_name = params["variables"][variable_i]["variable_name"].GetString()
		if (variable_name == "KratosMultiphysics.IGAStructuralMechanicsApplication.LOAD_TYPE"):
			condition_type = params["variables"][variable_i]["variable"].GetString()
			if (condition_type == "EDGE_LOAD"):
				conditionTypeInt = 1
			if (condition_type == "SURFACE_DEAD"):
				conditionTypeInt = 10
			if (condition_type == "SURFACE_PRESSURE"):
				conditionTypeInt = 100
			variables.update({variable_name : conditionTypeInt})
		if (variable_name == "KratosMultiphysics.IGAStructuralMechanicsApplication.DISTRIBUTED_LOAD_FACTOR"):
			variables.update({variable_name : params["variables"][variable_i]["variable"].GetDouble()})
		if (variable_name == "KratosMultiphysics.DIRECTION"):
			direction = KratosMultiphysics.Vector(3)
			direction[0] = params["variables"][variable_i]["variable"]["x"].GetDouble()
			direction[1] = params["variables"][variable_i]["variable"]["y"].GetDouble()
			direction[2] = params["variables"][variable_i]["variable"]["z"].GetDouble()
			variables.update({variable_name : direction})


	return IGAApplyLoad(model_part, variables)