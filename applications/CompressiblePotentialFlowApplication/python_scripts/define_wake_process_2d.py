import KratosMultiphysics
import KratosMultiphysics.CompressiblePotentialFlowApplication as CPFApp

def Factory(settings, Model):
    if(not isinstance(settings, KratosMultiphysics.Parameters)):
        raise Exception(
            "expected input shall be a Parameters object, encapsulating a json string")

    return DefineWakeProcess2D(Model, settings["Parameters"])
class DefineWakeProcess2D(KratosMultiphysics.Process):
    def __init__(self, Model, settings):
        # Call the base Kratos process constructor
        KratosMultiphysics.Process.__init__(self)

        # Check default settings
        default_settings = KratosMultiphysics.Parameters(r'''{
            "model_part_name": "",
            "epsilon": 1e-9,
            "compute_wake_at_each_step": false,
            "echo_level": 1
        }''')
        settings.ValidateAndAssignDefaults(default_settings)

        body_model_part_name = settings["model_part_name"].GetString()
        if body_model_part_name == "":
            err_msg = "Empty model_part_name in DefineWakeProcess2D\n"
            err_msg += "Please specify the model part that contains the body surface nodes"
            raise Exception(err_msg)
        self.body_model_part = Model[body_model_part_name]

        self.epsilon = settings["epsilon"].GetDouble()
        self.echo_level = settings["echo_level"].GetInt()

        self.fluid_model_part = self.body_model_part.GetRootModelPart()

        self.compute_wake_at_each_step = settings["compute_wake_at_each_step"].GetBool()

        for cond in self.body_model_part.Conditions:
            for node in cond.GetNodes():
                node.Set(KratosMultiphysics.SOLID)

    def ExecuteInitialize(self):

        CPFApp.Define2DWakeProcess(self.body_model_part, self.epsilon).ExecuteInitialize()

    def ExecuteInitializeSolutionStep(self):

        if self.compute_wake_at_each_step and self.fluid_model_part.ProcessInfo[KratosMultiphysics.STEP] > 1:
            CPFApp.Define2DWakeProcess(self.body_model_part, self.epsilon).ExecuteInitialize()

    def ExecuteFinalizeSolutionStep(self):
        if not self.fluid_model_part.HasSubModelPart("wake_sub_model_part"):
            raise Exception("Fluid model part does not have a wake_sub_model_part")
        else: self.wake_sub_model_part = self.fluid_model_part.GetSubModelPart("wake_sub_model_part")

        absolute_tolerance = 1e-9
        CPFApp.PotentialFlowUtilities.CheckIfWakeConditionsAreFulfilled2D(self.wake_sub_model_part, absolute_tolerance, self.echo_level)
        CPFApp.PotentialFlowUtilities.ComputePotentialJump2D(self.wake_sub_model_part)
