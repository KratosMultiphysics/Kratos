import KratosMultiphysics
import KratosMultiphysics.PoromechanicsApplication as KratosPoro

def Factory(settings, Model):
    if(type(settings) != KratosMultiphysics.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return LagrangianRotationProcess(Model, settings["Parameters"])

## All the processes python should be derived from "Process"

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

        # from math import pi
        # varying_parameters["angular_velocity_magnitude"] = - 2 * pi
        # varying_parameters["frame_rotation_axis_initial_point"] = [0., 0., 0.]
        # varying_parameters["frame_rotation_axis_final_point"] = [0., 0., 1.]

        self.process.append(KratosMultiphysics.ApplyConstantScalarValueProcess(model_part, params))

    def ExecuteInitialize(self):

        self.process.ExecuteInitialize()

    def ExecuteInitializeSolutionStep(self):

        self.process.ExecuteInitializeSolutionStep()
