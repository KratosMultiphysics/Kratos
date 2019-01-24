import KratosMultiphysics
import KratosMultiphysics.ShallowWaterApplication as KratosShallow

def Factory(settings, Model):
    if(type(settings) != KratosMultiphysics.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return ApplySlipProcess(Model, settings["Parameters"])

## All the processes python should be derived from "Process"
class ApplySlipProcess(KratosMultiphysics.Process):
    def __init__(self, Model, settings ):
        KratosMultiphysics.Process.__init__(self)

        default_parameters = KratosMultiphysics.Parameters("""
            {
                "model_part_name":"PLEASE_CHOOSE_MODEL_PART_NAME",
                "mesh_id": 0
            }
            """
            )
        settings.ValidateAndAssignDefaults(default_parameters)

        self.model_part = Model[settings["model_part_name"].GetString()]
        self.domain_size = self.model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE]

    def ExecuteInitialize(self):
        # Compute the normal on the nodes of interest
        KratosMultiphysics.NormalCalculationUtils().CalculateOnSimplex(self.model_part, self.domain_size)

        #TODO: Remove the IS_STRUCTURE variable set as soon as the flag SLIP migration is done
        #TODO: Remove the MESH_VELOCITY variable as soon the scheme doesn't use it
        for node in self.model_part.Nodes:
            node.SetValue(KratosMultiphysics.IS_STRUCTURE,1.0)
            node.SetSolutionStepValue(KratosMultiphysics.IS_STRUCTURE,0,1.0)
            node.SetSolutionStepValue(KratosMultiphysics.MESH_VELOCITY,0,[0.0, 0.0, 0.0])
