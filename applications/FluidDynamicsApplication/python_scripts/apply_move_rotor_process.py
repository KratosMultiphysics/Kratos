import KratosMultiphysics
import KratosMultiphysics.FluidDynamicsApplication as KratosFluid
import math

def Factory(settings, Model):
    if(type(settings) != KratosMultiphysics.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return ApplyMoveRotorProcess(Model, settings["Parameters"])

class ApplyMoveRotorProcess(KratosMultiphysics.Process):

    def __init__(self, Model, settings):

        KratosMultiphysics.Process.__init__(self)

        default_parameters = KratosMultiphysics.Parameters( """
        {
            "model_part_name"                                   :"rotor",
            "angular_velocity_with_respect_to_stator_center"    :0.0,
            "coordinates_of_stator_center_x"                    :0.0,
            "coordinates_of_stator_center_y"                    :0.0,
            "initial_coordinates_of_rotor_center_x"             :0.0,
            "initial_coordinates_of_rotor_center_y"             :0.0,
            "number_of_rotor_lobules"                            :0
        }  """ )

        settings.ValidateAndAssignDefaults(default_parameters);

        self.rotor_model_part = Model[settings["model_part_name"].GetString()]
        self.MoveRotorProcess = KratosFluid.MoveRotorProcess(self.rotor_model_part, settings)


    def ExecuteInitializeSolutionStep(self):
        self.MoveRotorProcess.ExecuteInitializeSolutionStep()
