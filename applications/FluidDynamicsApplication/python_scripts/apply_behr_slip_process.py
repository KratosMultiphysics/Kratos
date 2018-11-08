import KratosMultiphysics
import KratosMultiphysics.FluidDynamicsApplication

def Factory(settings, Model):

    if(type(settings) != KratosMultiphysics.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return BehrWallCondition(Model, settings["Parameters"])

## All the processes python should be derived from "Process"
class BehrWallCondition(KratosMultiphysics.Process):
    def __init__(self, Model, settings ):

        KratosMultiphysics.Process.__init__(self)

        ## Same as for regular slip
        default_parameters = KratosMultiphysics.Parameters( """
            {
                "model_part_name":"PLEASE_CHOOSE_MODEL_PART_NAME",
                "mesh_id": 0,
                "avoid_recomputing_normals": false
            }  """ );

        settings.ValidateAndAssignDefaults(default_parameters);

        self.model_part = Model[settings["model_part_name"].GetString()]
        self.avoid_recomputing_normals = settings["avoid_recomputing_normals"].GetBool()

        # Compute the normal on the nodes of interest -
        # Note that the model part employed here is supposed to only have slip "conditions"
        KratosMultiphysics.NormalCalculationUtils().CalculateOnSimplex(self.model_part, self.model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE])

        # Mark the nodes and conditions with the appropriate slip flag
        #TODO: Remove the IS_STRUCTURE variable set as soon as the flag SLIP migration is done
        for condition in self.model_part.Conditions: #TODO: this may well not be needed!
            
            condition.Set(KratosMultiphysics.SLIP, True)
            
            condition.SetValue(KratosMultiphysics.IS_STRUCTURE,1.0)
            # condition.Set(KratosMultiphysics.SLIP, False)

        #TODO: Remove the IS_STRUCTURE variable set as soon as the flag SLIP migration is done
        for node in self.model_part.Nodes:
            
            node.Set(KratosMultiphysics.SLIP, True)
            
            node.SetValue(KratosMultiphysics.IS_STRUCTURE,1.0)
            node.SetSolutionStepValue(KratosMultiphysics.IS_STRUCTURE,0,1.0)
            
            node.SetValue(KratosMultiphysics.Y_WALL,0.0)
            # node.Set(KratosMultiphysics.SLIP, False)


    def ExecuteInitializeSolutionStep(self):
        # Recompute the normals if needed (avoid if possible)
        if self.avoid_recomputing_normals == False:
            KratosMultiphysics.NormalCalculationUtils().CalculateOnSimplex(self.model_part, self.model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE])
            