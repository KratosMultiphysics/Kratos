import KratosMultiphysics
import KratosMultiphysics.ShallowWaterApplication as KratosShallow

def Factory(settings, Model):
    if(type(settings) != KratosMultiphysics.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return ApplySlipProcess(Model, settings["Parameters"])

##all the processes python processes should be derived from "python_process"
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

        # Compute the normal on the nodes of interest
        KratosMultiphysics.NormalCalculationUtils().CalculateOnSimplex(self.model_part, self.domain_size)

        # Mark the nodes and conditions with the appropriate slip flag
        #TODO: Remove the IS_STRUCTURE variable set as soon as the flag SLIP migration is done
        for condition in self.model_part.Conditions: #TODO: this may well not be needed!
            condition.SetValue(KratosMultiphysics.IS_STRUCTURE,1.0)

        #TODO: Remove the IS_STRUCTURE variable set as soon as the flag SLIP migration is done
        #TODO: Remove the MESH_VELOCITY variable as soon the scheme doesn't use it
        vel = [0.0, 0.0, 0.0]
        for node in self.model_part.Nodes:
            node.SetValue(KratosMultiphysics.IS_STRUCTURE,1.0)
            node.SetSolutionStepValue(KratosMultiphysics.IS_STRUCTURE,0,1.0)
            node.SetSolutionStepValue(KratosMultiphysics.MESH_VELOCITY,0,vel)
