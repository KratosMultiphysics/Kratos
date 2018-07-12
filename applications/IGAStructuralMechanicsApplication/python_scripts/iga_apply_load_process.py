import KratosMultiphysics


## All the processes python should be derived from "Process"
class IGAApplyLoad(KratosMultiphysics.Process):
    def __init__(self, model_part, variable_name, factor, direction, condition_type, mesh_id=0 ):
        KratosMultiphysics.Process.__init__(self) 
        
        variable = globals().get(variable_name)
        #print(model_part)
        for condition in model_part.GetMesh(mesh_id).Conditions:
            conditionTypeInt = 0;

            condition.SetValue(KratosMultiphysics.IGAStructuralMechanicsApplication.DISTRIBUTED_LOAD_FACTOR,factor)
            #condition.SetValue(KratosMultiphysics.IGAStructuralMechanicsApplication.CONDITION_TYPE_DEFINITION,condition_type)

            if (condition_type == "EDGE_LOAD"):
                conditionTypeInt = 1
                condition.SetValue(KratosMultiphysics.DIRECTION,[direction[0].GetDouble(),direction[1].GetDouble(),direction[2].GetDouble()])

            if (condition_type == "SURFACE_DEAD"):
                conditionTypeInt = 10
                condition.SetValue(KratosMultiphysics.DIRECTION,[direction[0].GetDouble(),direction[1].GetDouble(),direction[2].GetDouble()])

            if (condition_type == "SURFACE_PRESSURE"):
                conditionTypeInt = 100

            condition.SetValue(KratosMultiphysics.IGAStructuralMechanicsApplication.LOAD_TYPE,conditionTypeInt)

            #print("Finished construction of IGAApplyLoad Process")
        
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
	factor = params["factor"].GetDouble()
	direction = params["direction"]

	condition_type = settings["condition_type_description"].GetString()

	return IGAApplyLoad(model_part, variable_name, factor, direction, condition_type)