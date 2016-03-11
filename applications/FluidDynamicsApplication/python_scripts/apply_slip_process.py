import KratosMultiphysics 
import KratosMultiphysics.FluidDynamicsApplication 

def Factory(settings, Model):
    if(type(settings) != KratosMultiphysics.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return ApplySlipProcess(Model, settings["Parameters"])

##all the processes python processes should be derived from "python_process"
class ApplySlipProcess(KratosMultiphysics.Process):
    def __init__(self, Model, settings ):
        KratosMultiphysics.Process.__init__(self)
        
        default_parameters = KratosMultiphysics.Parameters( """
            {
                "model_part_name":"PLEASE_CHOOSE_MODEL_PART_NAME",
                "mesh_id": 0,
                "avoid_recomputing_normals": false
            }  """ );
        
        settings.ValidateAndAssignDefaults(default_parameters);
        
        self.model_part = Model[settings["model_part_name"].GetString()]
        self.avoid_recomputing_normals = settings["avoid_recomputing_normals"].GetBool()

        #compute the normal on the nodes of interest - note that the model part employed here is supposed to
        #only have slip "conditions"
        KratosMultiphysics.NormalCalculationUtils().CalculateOnSimplex(self.model_part, self.model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE])
        
        #mark the nodes and conditions with the appropriate slip flag
        #TODO: a flag shall be used here!!!!!
        for cond in self.model_part.Conditions: #TODO: this may well not be needed!
            cond.SetValue(KratosMultiphysics.IS_STRUCTURE,1.0)
            
        #TODO: use a flag!!!
        KratosMultiphysics.VariableUtils().SetScalarVar(KratosMultiphysics.IS_STRUCTURE,1.0, self.model_part.Nodes)
        
    def InitializeSolutionStep(self):
        #recompute the normals
        if self.avoid_recomputing_normals == False:
            NormalCalculationUtils().CalculateOnSimplex(self.model_part, self.model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE])
