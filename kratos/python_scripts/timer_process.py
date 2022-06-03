# Importing the Kratos Library
import KratosMultiphysics

def Factory(settings, Model):
    if not isinstance(settings, KratosMultiphysics.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return TimerProcess(Model, settings["Parameters"])

# All the processes python processes should be derived from "Process"
class TimerProcess(KratosMultiphysics.Process):
    """This process helps to measure the time consumed on the simulations

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
            "help"                       : "This process helps to measure the time consumed on the simulations",
            "output_filename"            : "",
            "print_interval_information" : false,
            "interval_name"              : "Analysis"
        }
        """
        )

        # Assign this here since it will change the "interval" prior to validation
        settings.ValidateAndAssignDefaults(default_settings)

        self.interval_name = settings["interval_name"].GetString()
        self.output_filename = settings["output_filename"].GetString()

        # Defining timer
        self.timer = KratosMultiphysics.Timer()

        # Interval information
        self.timer.SetPrintIntervalInformation(settings["print_interval_information"].GetBool())

        # Output file
        if self.output_filename != "":
            self.timer.SetOuputFile(self.output_filename)
        else:
            self.timer.SetPrintOnScreen(True)

        # Starting timer
        self.timer.Start(self.interval_name)

    def ExecuteFinalize(self):
        """ This function is designed for being called at the end of the computations

        Keyword arguments:
        self -- It signifies an instance of a class.
        """
        self.timer.Stop(self.interval_name)
        if self.output_filename != "":
            self.timer.PrintTimingInformation(self.timer)
            self.timer.CloseOuputFile()
