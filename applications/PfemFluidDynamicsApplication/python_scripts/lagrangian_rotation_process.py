import KratosMultiphysics
import KratosMultiphysics.PfemFluidDynamicsApplication as KratosPfemFluid

def Factory(settings, Model):
    if(type(settings) != KratosMultiphysics.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return LagrangianRotationProcess(Model, settings["Parameters"])

class LagrangianRotationProcess(KratosMultiphysics.Process):
    def __init__(self, Model, settings ):
        KratosMultiphysics.Process.__init__(self)

        model_part = Model[settings["model_part_name"].GetString()]

        params = KratosMultiphysics.Parameters("{}")
        params.AddValue("model_part_name",settings["model_part_name"])
        params.AddValue("angular_velocity",settings["angular_velocity"])
        params.AddValue("rotation_axis_initial_point",settings["rotation_axis_initial_point"])
        params.AddValue("rotation_axis_final_point",settings["rotation_axis_final_point"])
        params.AddValue("initial_time",settings["initial_time"])

        self.process = KratosPfemFluid.LagrangianRotationProcess(model_part, params)

    def ExecuteInitialize(self):

        self.process.ExecuteInitialize()

    def ExecuteInitializeSolutionStep(self):

        self.process.ExecuteInitializeSolutionStep()
