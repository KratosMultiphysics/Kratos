import KratosMultiphysics
import KratosMultiphysics.SwimmingDEMApplication as KratosSDEM

def Factory(settings, Model):
    if(type(settings) != KratosMultiphysics.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return ApplyRigidRotationProcess(Model, settings["Parameters"])

class ApplyRigidRotationProcess(KratosMultiphysics.Process):
    def __init__(self, Model, settings ):
        KratosMultiphysics.Process.__init__(self)

        # Fluid walls model part
        self.model_part = Model[settings["model_part_name"].GetString()]

        self.params = KratosMultiphysics.Parameters("{}")
        self.params.AddValue("model_part_name",settings["model_part_name"])
        self.params.AddValue("angular_velocity",settings["angular_velocity"])
        self.params.AddValue("rotation_axis_initial_point",settings["rotation_axis_initial_point"])
        self.params.AddValue("rotation_axis_final_point",settings["rotation_axis_final_point"])
        self.params.AddValue("initial_time",settings["initial_time"])

        self.process = KratosSDEM.ApplyRigidRotationProcess(self.model_part, self.params)

    def ExecuteInitialize(self):

        self.process.ExecuteInitialize()

    def ExecuteInitializeSolutionStep(self):

        self.process.ExecuteInitializeSolutionStep()
