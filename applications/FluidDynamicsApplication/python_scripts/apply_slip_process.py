import KratosMultiphysics
import KratosMultiphysics.FluidDynamicsApplication

def Factory(settings, Model):
    if(type(settings) != KratosMultiphysics.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return ApplySlipProcess(Model, settings["Parameters"])

## All the processes python should be derived from "Process"
class ApplySlipProcess(KratosMultiphysics.Process):
    def __init__(self, Model, settings ):
        KratosMultiphysics.Process.__init__(self)

        default_parameters = KratosMultiphysics.Parameters( """
            {
                "model_part_name":"PLEASE_CHOOSE_MODEL_PART_NAME",
                "mesh_id": 0,
                "avoid_recomputing_normals": false,
                "uniform_navier_slip_length" : 0.01
            }  """ )

        self.navier_slip_active = False
        if settings.Has("uniform_navier_slip_length"):
            self.navier_slip_active = True

        settings.ValidateAndAssignDefaults(default_parameters)

        self.model_part = Model[settings["model_part_name"].GetString()]
        self.avoid_recomputing_normals = settings["avoid_recomputing_normals"].GetBool()

        # Mark the nodes and conditions with the appropriate slip flag
        for condition in self.model_part.Conditions:
            condition.Set(KratosMultiphysics.SLIP, True)

        for node in self.model_part.Nodes:
            node.Set(KratosMultiphysics.SLIP, True)
            node.SetValue(KratosMultiphysics.Y_WALL,0.0)

        if self.navier_slip_active:
            navier_slip_length = settings["uniform_navier_slip_length"].GetDouble()
            KratosMultiphysics.VariableUtils().SetNonHistoricalVariable(
                KratosMultiphysics.FluidDynamicsApplication.SLIP_LENGTH,
                navier_slip_length,
                self.model_part.Nodes)
        else:
            KratosMultiphysics.VariableUtils().SetNonHistoricalVariable(
                KratosMultiphysics.FluidDynamicsApplication.SLIP_LENGTH,
                1.0e8,
                self.model_part.Nodes)

    def ExecuteInitialize(self):
        # Compute the normal on the nodes of interest -
        # Note that the model part employed here is supposed to only have slip "conditions"
        enforce_generic_algorithm = True
        KratosMultiphysics.NormalCalculationUtils().CalculateNormals(self.model_part, enforce_generic_algorithm)

    def ExecuteInitializeSolutionStep(self):
        # Recompute the normals if needed
        if self.avoid_recomputing_normals == False:
            enforce_generic_algorithm = True
            KratosMultiphysics.NormalCalculationUtils().CalculateNormals(self.model_part, enforce_generic_algorithm)
