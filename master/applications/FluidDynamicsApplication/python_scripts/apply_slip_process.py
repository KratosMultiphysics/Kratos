import KratosMultiphysics
import KratosMultiphysics.FluidDynamicsApplication
from KratosMultiphysics.kratos_utilities import IssueDeprecationWarning

def Factory(settings, Model):
    if(type(settings) != KratosMultiphysics.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return ApplySlipProcess(Model, settings["Parameters"])

## All the processes python should be derived from "Process"
class ApplySlipProcess(KratosMultiphysics.Process):
    def __init__(self, Model, settings ):
        KratosMultiphysics.Process.__init__(self)

        # TODO: Wipe deprecated options. To be removed after the deprecation period.
        if settings.Has("uniform_navier_slip_length"):
            settings.RemoveValue("uniform_navier_slip_length")
            IssueDeprecationWarning('ApplySlipProcess', 'ApplySlipProcess no longer supports the Navier-slip wall behavior. Check the ApplyFluidWallProcess.')

        # Validate and assign settings
        settings.ValidateAndAssignDefaults(self.GetDefaultParameters())
        if not settings["model_part_name"].GetString():
            raise Exception("'model_part_name' is not provided.")
        self.model_part = Model[settings["model_part_name"].GetString()]
        self.avoid_recomputing_normals = settings["avoid_recomputing_normals"].GetBool()
        self.model_part.ProcessInfo[KratosMultiphysics.FluidDynamicsApplication.SLIP_TANGENTIAL_CORRECTION_SWITCH] = settings["slip_tangential_correction"].GetBool()

    @classmethod
    def GetDefaultParameters(cls):
        default_parameters = KratosMultiphysics.Parameters("""{
            "model_part_name" : "",
            "avoid_recomputing_normals" : false,
            "slip_tangential_correction" : true
        }""")

        return default_parameters

    def ExecuteInitialize(self):
        # Mark the nodes and conditions with the SLIP flag
        KratosMultiphysics.VariableUtils().SetFlag(KratosMultiphysics.SLIP, True, self.model_part.Nodes)
        KratosMultiphysics.VariableUtils().SetFlag(KratosMultiphysics.SLIP, True, self.model_part.Conditions)

        # Compute the normal on the nodes of interest -
        # Note that the model part employed here is supposed to only have slip "conditions"
        enforce_generic_algorithm = True
        KratosMultiphysics.NormalCalculationUtils().CalculateNormals(self.model_part, enforce_generic_algorithm)

    def ExecuteInitializeSolutionStep(self):
        # Recompute the normals if needed
        if self.avoid_recomputing_normals == False:
            enforce_generic_algorithm = True
            KratosMultiphysics.NormalCalculationUtils().CalculateNormals(self.model_part, enforce_generic_algorithm)
