# Importing the Kratos Library
import KratosMultiphysics as KM

# Define the process factory
def Factory(settings, Model):
    if not isinstance(settings, KM.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return MoveMeshProcess(Model, settings["Parameters"])

class MoveMeshProcess(KM.Process):
    """
    This process moves the mesh of a model part using the specified variable.
    The mesh is moved only if the current time is within the specified interval.

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
        # call the base class constructor
        super().__init__()

        # Validate input settings
        default_settings = KM.Parameters("""
        {
            "help"            : "This process moves the mesh of a model part using the specified variable. The mesh is moved only if the current time is within the specified interval.",
            "model_part_name" : "please_specify_model_part_name",
            "variable_name"   : "DISPLACEMENT",
            "interval"        : [0.0, 1e30]
        }
        """
        )

        # Validate and assign
        settings.ValidateAndAssignDefaults(default_settings)

        # Define member variables of the class
        self.model_part = Model[settings["model_part_name"].GetString()]
        self.variable = KM.KratosGlobals.GetVariable(settings["variable_name"].GetString())
        self.interval = KM.IntervalUtility(settings)

    def ExecuteInitializeSolutionStep(self):
        """ This method is executed in order to initialize the current step

        Keyword arguments:
        self -- It signifies an instance of a class.
        """
        process_info = self.model_part.ProcessInfo
        current_time = process_info[KM.TIME]

        # Check if the current time is in the interval
        if self.interval.IsInInterval(current_time):
            # If the interval is valid, we apply the mesh moving process
            KM.BaseSolvingStrategy.SimpleDisplacementMeshMoving(self.model_part, self.variable)

