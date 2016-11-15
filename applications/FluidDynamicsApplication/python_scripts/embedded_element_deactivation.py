import KratosMultiphysics

def Factory(settings, Model):
    if(type(settings) != KratosMultiphysics.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return EmbeddedElementDeactivation(Model, settings["Parameters"])

class EmbeddedElementDeactivation(KratosMultiphysics.Process):
    
    def __init__(self, Model, settings):
        
        KratosMultiphysics.Process.__init__(self)

        default_parameters = KratosMultiphysics.Parameters( """
        {
            "mesh_id"                   : 0,
            "model_part_name"           : "CHOOSE_FLUID_MODELPART_NAME",
            "check_at_each_time_step"   : false
        }  """ )

        settings.ValidateAndAssignDefaults(default_parameters);
        
        self.fluid_model_part = Model[settings["model_part_name"].GetString()]
        self.nnodes = self.fluid_model_part.ProcessInfo.GetValue(KratosMultiphysics.DOMAIN_SIZE) + 1 # Note that this process only works for simplicial elements (tetrahedras and triangles)
        self.check_at_each_time_step = settings["check_at_each_time_step"].GetBool()
                
                
    def ExecuteInitialize(self):
        
        for element in self.fluid_model_part.Elements:
            inside = 0
            
            # Check the sign of the distance function at the element nodes
            for node in element.GetNodes():
                if(node.GetSolutionStepValue(KratosMultiphysics.DISTANCE) < 0.0):
                    inside += 1
                    
            # If all the nodes have negative distance value (structure domain) deactivate the element
            if(inside == self.nnodes):
                element.Set(KratosMultiphysics.ACTIVE,False)
                
                
    def ExecuteInitializeSolutionStep(self):
        
        if (self.check_at_each_time_step == True):            
                        
            for element in self.fluid_model_part.Elements:
                
                inside = 0
                
                # Check the sign of the distance function at the element nodes
                for node in element.GetNodes():
                    if(node.GetSolutionStepValue(KratosMultiphysics.DISTANCE) < 0.0):
                        inside += 1
                        
                # If all the nodes have negative distance value (structure domain) deactivate the element
                if(inside == self.nnodes):
                    element.Set(KratosMultiphysics.ACTIVE, False)
                # Otherwise, activate the element again
                elif(inside < self.nnodes):
                    element.Set(KratosMultiphysics.ACTIVE, True)
                # Security check
                else:
                    raise Exception ("ERROR: The number of inside elements is larger than the element nodes!")
        
        else:
            pass
