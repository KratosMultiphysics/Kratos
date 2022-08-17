# Importing the Kratos Library
import KratosMultiphysics
import KratosMultiphysics.StructuralMechanicsApplication as KSM

def Factory(settings, Model):
    if not isinstance(settings, KratosMultiphysics.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return SetMovingLoadProcess(Model, settings["Parameters"])

from KratosMultiphysics.StructuralMechanicsApplication import set_moving_load_process as smlp
#
## All the processes python should be derived from "Process"
class SetMovingLoadProcess(KratosMultiphysics.Process):
    """This process assigns given POINT_LOAD and PRESCRIBED_DISPLACEMENT values
    (scalar) to the DisplacementControl belonging a certain submodelpart.
    Currently, it is restricted to one node and one global direction.

    Only the member variables listed below should be accessed directly.

    Public member variables:
    Model -- the container of the different model parts.
    settings -- Kratos parameters containing solver settings.
    """

    def __init__(self, Model, settings ):
        """ The default constructor of the class

        Keyword arguments:
        self -- It signifies an instance of a class.
        Model -- the container of the different model parts.
        settings -- Kratos parameters containing solver settings.
        """

        KratosMultiphysics.Process.__init__(self)

        #The value can be a double or a string (function)
        default_settings = KratosMultiphysics.Parameters("""
        {
            "help"            : "This process applies a moving load condition belonging to a modelpart. The load moves over line elements.",
            "model_part_name" : "please_specify_model_part_name",
            "variable_name"   : "MOVING_LOAD",
            "load"            : [0.0, 1.0, 0.0],
            "direction"       : [1,1,1],
            "velocity"        : 1
        }
        """
        )

        settings.ValidateAndAssignDefaults(default_settings)


        params_load = KratosMultiphysics.Parameters("{}")
        params_load.AddValue("model_part_name", settings["model_part_name"])
        params_load.AddValue("variable_name", settings["variable_name"])
        params_load.AddValue("load", settings["load"])
        params_load.AddValue("direction", settings["direction"])
        params_load.AddValue("velocity", settings["velocity"])

        # Set process
        self.model_part = Model.GetModelPart(params_load["model_part_name"].GetString())
        self.process = KSM.SetMovingLoadProcess(self.model_part, params_load)

    def ExecuteInitialize(self):
        """ This method is executed in order to initialize the current step

        Keyword arguments:
        self -- It signifies an instance of a class.
        """

        self.process.ExecuteInitialize()

    def ExecuteInitializeSolutionStep(self):
    #     """ This method is executed in order to initialize the current step
    #
    #     Keyword arguments:
    #     self -- It signifies an instance of a class.
    #     """
        self.process.ExecuteInitializeSolutionStep()

    def ExecuteFinalizeSolutionStep(self):
    #     """ This method is executed in order to initialize the current step
    #
    #     Keyword arguments:
    #     self -- It signifies an instance of a class.
    #     """
        self.process.ExecuteFinalizeSolutionStep()