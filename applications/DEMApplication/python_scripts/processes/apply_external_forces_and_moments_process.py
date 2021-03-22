# Importing the Kratos Library
import KratosMultiphysics
import KratosMultiphysics.DEMApplication as DEM

def Factory(settings, Model):
    if not isinstance(settings, KratosMultiphysics.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return ApplyExternalForcesAndMomentsProcess(Model, settings["Parameters"])

## All the processes python should be derived from "Process"
class ApplyExternalForcesAndMomentsProcess(KratosMultiphysics.Process):
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
            "force_settings" : {
                "value"                : [10.0, "3*t", "x+y"]
            },
            "moment_settings" : {
                "value"                : [10.0, "3*t", "x+y"]
            },
            "interval"             : [0.0, 1e30]
        }
        """
        )
        #example of admissible values for "value" : [10.0, "3*t", "x+y"]

        settings.ValidateAndAssignDefaults(default_settings)

        self.model_part = Model[settings["model_part_name"].GetString()]
        self.cplusplus_version_of_this_process = DEM.ApplyExternalForcesAndMomentsProcess(self.model_part, settings)


    def ExecuteInitializeSolutionStep(self):
        self.cplusplus_version_of_this_process.ExecuteInitializeSolutionStep()

    def ExecuteFinalizeSolutionStep(self):
        self.cplusplus_version_of_this_process.ExecuteFinalizeSolutionStep()

