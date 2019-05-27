import KratosMultiphysics
import KratosMultiphysics.DrillingBeamApplication as DBA
import math

def Factory(settings, Model):
    if(type(settings) != KratosMultiphysics.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return ApplyAdvanceAndRotatePartOfStructureProcess(Model, settings["Parameters"])

class ApplyAdvanceAndRotatePartOfStructureProcess(KratosMultiphysics.Process):

    def __init__(self, Model, settings):

        KratosMultiphysics.Process.__init__(self)

        default_parameters = KratosMultiphysics.Parameters( """
        {
                "model_part_name":"to_be_rotated",
                "advance_velocity": 0.0,
                "angular_velocity_x":0.0,
                "angular_velocity_y":0.0,
                "angular_velocity_z":0.0,
                "coordinates_rotation_center_x":0.0,
                "coordinates_rotation_center_y":0.0,
                "coordinates_rotation_center_z":0.0
        }  """ )

        settings.ValidateAndAssignDefaults(default_parameters);

        self.model_part_to_be_rotated = Model[settings["model_part_name"].GetString()]
        self.AdvanceAndRotatePartOfStructureProcess = DBA.AdvanceAndRotatePartOfStructureProcess(self.model_part_to_be_rotated, settings)


    def ExecuteInitializeSolutionStep(self):
        self.AdvanceAndRotatePartOfStructureProcess.ExecuteInitializeSolutionStep()
