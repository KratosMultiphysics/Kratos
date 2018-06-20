import KratosMultiphysics


## All the processes python should be derived from "Process"
class IGAApplySupport(KratosMultiphysics.Process):
    def __init__(self, model_part, variable_name, is_fixed_x, is_fixed_y, is_fixed_z, is_fixed_rot, value, mesh_id=0 ):
        KratosMultiphysics.Process.__init__(self) 
       
        DisplacementRotationFix = 0 #defined by rot, dispx, dispy, dispz
        if (is_fixed_rot):
            DisplacementRotationFix += 1000
        if (is_fixed_x):
            DisplacementRotationFix += 100
        if (is_fixed_y):
            DisplacementRotationFix += 10
        if (is_fixed_z):
            DisplacementRotationFix += 1

        variable = globals().get(variable_name)
        #print(model_part)
        for condition in model_part.GetMesh(mesh_id).Conditions:
            condition.SetValue(KratosMultiphysics.DISPLACEMENT,[value[0].GetDouble(),value[1].GetDouble(),value[2].GetDouble()])
            condition.SetValue(KratosMultiphysics.IGAStructuralMechanicsApplication.DISPLACEMENT_ROTATION_FIX, DisplacementRotationFix)
        #print(value[0].GetDouble(),value[1].GetDouble(),value[2].GetDouble())
        print("Finished construction of IGAApplySupport Process")
        
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
	model_part = Model.get(  params["model_part_name"].GetString() , "model part not found" )
	mesh_id = params["mesh_id"]

	#if(settings["process_name"] == "IGAApplyLoad"):
	variable_name = params["variable_name"] 
	is_fixed_x = params["is_fixed_x"].GetBool()
	is_fixed_y = params["is_fixed_y"].GetBool()
	is_fixed_z = params["is_fixed_z"].GetBool()
	is_fixed_rot = params["is_fixed_rot"].GetBool()
	value = params["value"]
	
	return IGAApplySupport(model_part, variable_name, is_fixed_x, is_fixed_y, is_fixed_z, is_fixed_rot, value)