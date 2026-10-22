import KratosMultiphysics as KM
import KratosMultiphysics.FluidDynamicsApplication as KM_CFD

from KratosMultiphysics.assign_vector_variable_process import AssignVectorVariableProcess


def Factory(settings, Model):
    if(type(settings) != KM.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return ApplyNitscheNavierSlipProcess(Model, settings["Parameters"])


## All the processes python should be derived from "Process"
class ApplyNitscheNavierSlipProcess(KM.Process):
    def __init__(self, Model, settings ):
        KM.Process.__init__(self)

        # Validate and assign settings
        settings.ValidateAndAssignDefaults(self.GetDefaultParameters())
        if not settings["model_part_name"].GetString():
            raise Exception("'model_part_name' is not provided.")
        self.model_part = Model[settings["model_part_name"].GetString()]
        self.avoid_recomputing_normals = settings["avoid_recomputing_normals"].GetBool()

        self.slip_length = settings["slip_length"].GetDouble()
        self.penalty_coefficient = settings["penalty_coefficient"].GetDouble()

    @classmethod
    def GetDefaultParameters(cls):
        default_parameters = KM.Parameters("""{
            "model_part_name" : "",
            "avoid_recomputing_normals" : false,
            "slip_length" : 0.0,
            "penalty_coefficient" : 10.0,
            "variable_name"        : "VELOCITY",
            "interval"             : [0.0, 1e30],
            "value"                : [0.0, 0.0, 0.0],
            "constrained"          : [true,true,true],
            "local_axes"           : {}
        }""")

        return default_parameters

    def ExecuteInitialize(self):
        # Get the condition registry name
        # NOTE that in here we are assuming a unique condition geometry in model part
        pts_num = 0
        for condition in self.model_part.Conditions:
            pts_num = condition.GetGeometry().PointsNumber()
            break
        pts_num = self.model_part.GetCommunicator().GetDataCommunicator().MaxAll(pts_num)
        domain_size = self.model_part.ProcessInfo[KM.DOMAIN_SIZE]
        cond_reg_name = f'{"NavierStokesNitscheNavierSlipWallCondition"}{domain_size}D{pts_num}N'

        # Substitute current conditions by the corresponding ones implementing the wall model
        replace_settings = KM.Parameters("""{
            "element_name" : "",
            "condition_name" : ""
        }""")
        replace_settings["condition_name"].SetString(cond_reg_name)
        KM.ReplaceElementsAndConditionsProcess(self.model_part, replace_settings).Execute()

        # Compute the normal on the nodes of interest
        KM.NormalCalculationUtils().CalculateOnSimplex(self.model_part, domain_size)

        # Set the WALL flag in the wall model part new conditions
        # This is required in order to add the wall model contribution within the condition implementation
        KM.VariableUtils().SetFlag(KM.WALL, True, self.model_part.Conditions)
        KM.VariableUtils().SetFlag(KM.SLIP, False, self.model_part.Nodes)  #TODO not necessary?

        # Assign provided slip length to the nodes
        for condition in self.model_part.Conditions:
            for node in condition.GetNodes():
                node.SetValue(KM_CFD.SLIP_LENGTH, self.slip_length)
                node.SetValue(KM.PENALTY_COEFFICIENT, self.penalty_coefficient)
        self.model_part.GetCommunicator().SynchronizeNonHistoricalVariable(KM_CFD.SLIP_LENGTH)
        self.model_part.GetCommunicator().SynchronizeNonHistoricalVariable(KM.PENALTY_COEFFICIENT)

    def ExecuteInitializeSolutionStep(self):
        # Recompute the normals if needed (e.g. moving boundaries)
        if self.avoid_recomputing_normals == False:
            domain_size = self.model_part.ProcessInfo[KM.DOMAIN_SIZE]
            KM.NormalCalculationUtils().CalculateOnSimplex(self.model_part, domain_size)
