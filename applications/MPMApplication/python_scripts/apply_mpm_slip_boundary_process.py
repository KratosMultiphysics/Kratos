import KratosMultiphysics
import KratosMultiphysics.MPMApplication as KratosMPM

def Factory(settings, Model):
    if(not isinstance(settings, KratosMultiphysics.Parameters)):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return ApplyMPMSlipBoundaryProcess(Model, settings["Parameters"])

## All the processes python should be derived from "Process"
class ApplyMPMSlipBoundaryProcess(KratosMultiphysics.Process):
    def __init__(self, Model, settings ):
        KratosMultiphysics.Process.__init__(self)

        default_parameters = KratosMultiphysics.Parameters( """
            {
                "model_part_name"           :"PLEASE_CHOOSE_MODEL_PART_NAME",
                "mesh_id"                   : 0,
                "friction_coefficient"      : 0,
                "tangential_penalty_factor" : 0,
                "option"                    : "",
                "avoid_recomputing_normals" : true
            }  """ )

        settings.ValidateAndAssignDefaults(default_parameters)

        self.model_part = Model[settings["model_part_name"].GetString()]
        self.avoid_recomputing_normals = settings["avoid_recomputing_normals"].GetBool()

        # get friction parameters
        self.friction_coefficient = settings["friction_coefficient"].GetDouble()
        self.tangential_penalty_factor = settings["tangential_penalty_factor"].GetDouble()
        self.option = settings["option"].GetString()

        if self.friction_coefficient > 0 and self.tangential_penalty_factor > 0:
            # friction active -- set flag on ProcessInfo
            self.model_part.ProcessInfo[KratosMPM.FRICTION_ACTIVE] = True


        # Compute the normal on the nodes of interest -
        # Note that the model part employed here is supposed to only have slip "conditions"
        KratosMultiphysics.NormalCalculationUtils().CalculateOnSimplex(self.model_part, self.model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE])

        self.flip_normal = (self.option == "flip_normal")

        #TODO: Remove the IS_STRUCTURE variable set as soon as the flag SLIP migration is done
        for node in self.model_part.Nodes:
            node.Set(KratosMultiphysics.SLIP, True)
            node.SetValue(KratosMultiphysics.IS_STRUCTURE,1.0)
            node.SetSolutionStepValue(KratosMultiphysics.IS_STRUCTURE,0,1.0)
            node.SetValue(KratosMultiphysics.FRICTION_COEFFICIENT, self.friction_coefficient)
            node.SetValue(KratosMPM.TANGENTIAL_PENALTY_FACTOR, self.tangential_penalty_factor)
            node.Set(KratosMultiphysics.MODIFIED, self.flip_normal)



    def ExecuteInitializeSolutionStep(self):
        # Recompute the normals if needed
        if self.avoid_recomputing_normals == False:
            KratosMultiphysics.NormalCalculationUtils().CalculateOnSimplex(self.model_part, self.model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE])
