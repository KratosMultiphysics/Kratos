from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
# Importing the Kratos Library
import KratosMultiphysics
import KratosMultiphysics.DEMApplication as DEM

def Factory(settings, Model):
    if not isinstance(settings, KratosMultiphysics.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return ApplyKinematicConstraintsProcess(Model, settings["Parameters"])

## All the processes python should be derived from "Process"
class ApplyKinematicConstraintsProcess(KratosMultiphysics.Process):
    """This process assigns a given value (vector) to the nodes belonging a certain submodelpart

    Only the member variables listed below should be accessed directly.

    Public member variables:
    Model -- the container of the different model parts.
    settings -- Kratos parameters containing solver settings.
    """

    def __init__(self, Model, settings):
        """ The default constructor of the class

        Keyword arguments:
        self -- It signifies an instance of a class.
        Model -- the container of the different model parts.
        settings -- Kratos parameters containing solver settings.
        """

        KratosMultiphysics.Process.__init__(self)

        default_settings = KratosMultiphysics.Parameters("""
        {
            "help"                 : "This process applies constraints to the particles in a certain submodelpart, for a certain time interval",
            "mesh_id"              : 0,
            "model_part_name"      : "please_specify_model_part_name",
            "velocity_constraints_settings" : {
                "constrained"          : [true,true,true],
                "value"                : [10.0, "3*t", "x+y"]
            },
            "angular_velocity_constraints_settings" : {
                "constrained"          : [true,true,true],
                "value"                : [10.0, "3*t", "x+y"]
            },
            "interval"             : [0.0, 1e30]
        }
        """
        )
        #example of admissible values for "value" : [10.0, "3*t", "x+y"]

        ## Trick to ensure that if someone sets constrained as a single bool, it is transformed to a vector
        if settings["velocity_constraints_settings"].Has("constrained"):
            if settings["velocity_constraints_settings"]["constrained"].IsBool():
                is_fixed = settings["velocity_constraints_settings"]["constrained"].GetBool()
                #print("is_fixed = ",is_fixed)
                settings["velocity_constraints_settings"]["constrained"] = default_settings["velocity_constraints_settings"]["constrained"]
                for i in range(3):
                    settings["velocity_constraints_settings"]["constrained"][i].SetBool(is_fixed)

        if settings["angular_velocity_constraints_settings"].Has("constrained"):
            if settings["angular_velocity_constraints_settings"]["constrained"].IsBool():
                is_fixed = settings["angular_velocity_constraints_settings"]["constrained"].GetBool()
                #print("is_fixed = ",is_fixed)
                settings["angular_velocity_constraints_settings"]["constrained"] = default_settings["angular_velocity_constraints_settings"]["constrained"]
                for i in range(3):
                    settings["angular_velocity_constraints_settings"]["constrained"][i].SetBool(is_fixed)

        settings.ValidateAndAssignDefaults(default_settings)

        self.model_part = Model[settings["model_part_name"].GetString()]
        self.cplusplus_version_of_this_process = DEM.ApplyKinematicConstraintsProcess(self.model_part, settings)


    def ExecuteInitializeSolutionStep(self):
        self.cplusplus_version_of_this_process.ExecuteInitializeSolutionStep()

    def ExecuteFinalizeSolutionStep(self):
        self.cplusplus_version_of_this_process.ExecuteFinalizeSolutionStep()

