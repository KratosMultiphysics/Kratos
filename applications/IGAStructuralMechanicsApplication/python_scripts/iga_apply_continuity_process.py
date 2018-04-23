import KratosMultiphysics
import python_process

##all the processes python processes should be derived from "python_process"
class IGAApplyContinuity(python_process.PythonProcess):
    def __init__(self, model_part, variables, mesh_id=0 ):
        python_process.PythonProcess.__init__(self) 

        for condition in model_part.Conditions:
            for variable_key in variables:
                condition.SetValue(eval(variable_key), variables[variable_key])

        print("Finished construction of IGAApplyContinuity Process")
        
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
	model_part = Model[params["model_part_name"].GetString()].GetSubModelPart(params["sub_model_part_name"].GetString())# , "model part not found" )
	mesh_id = params["mesh_id"]

	variables = {}
	#if(settings["process_name"] == "IGAApplyLoad"):
	for variable_i in range (0,params["variables"].size()):
		variable_name = params["variables"][variable_i]["variable_name"].GetString()
		if (variable_name == "KratosMultiphysics.IGAStructuralMechanicsApplication.PENALTY_FACTOR"):
			variables.update({variable_name : params["variables"][variable_i]["variable"].GetDouble()})
		if (variable_name == "KratosMultiphysics.IGAStructuralMechanicsApplication.DISPLACEMENT_ROTATION_FIX"):
			DisplacementRotationFix = 0 #defined by rot, dispx, dispy, dispz
			if (params["variables"][variable_i]["variable"]["C1-Continuity"]["t1"].GetBool()):
				DisplacementRotationFix += 1000
			if (params["variables"][variable_i]["variable"]["C0-Continuity"]["x"].GetBool()):
				DisplacementRotationFix += 100
			if (params["variables"][variable_i]["variable"]["C0-Continuity"]["y"].GetBool()):
				DisplacementRotationFix += 10
			if (params["variables"][variable_i]["variable"]["C0-Continuity"]["z"].GetBool()):
				DisplacementRotationFix += 1
			variables.update({variable_name : DisplacementRotationFix})
		if (variable_name == "KratosMultiphysics.DISPLACEMENT"):
			displacements = KratosMultiphysics.Vector(3)
			displacements[0] = params["variables"][variable_i]["variable"]["C0-Continuity"]["x"].GetDouble()
			displacements[1] = params["variables"][variable_i]["variable"]["C0-Continuity"]["y"].GetDouble()
			displacements[2] = params["variables"][variable_i]["variable"]["C0-Continuity"]["z"].GetDouble()
			variables.update({variable_name : displacements})

	
	return IGAApplyContinuity(model_part, variables)