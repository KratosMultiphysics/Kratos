# Importing the Kratos Library
import KratosMultiphysics

# Import KratosUnittest
import KratosMultiphysics.KratosUnittest as KratosUnittest

def Factory(settings, model):
    if not isinstance(settings, KratosMultiphysics.Parameters):
        raise Exception("Expected input shall be a Parameters object, encapsulating a json string")
    return FromJsonCheckResultProcess(model, settings["Parameters"])

class FromJsonCheckResultProcess(KratosMultiphysics.Process, KratosUnittest.TestCase):
    """This class is used in order to check results using a json file
    containing the solution a given model part with a certain frequency

    Only the member variables listed below should be accessed directly.

    Public member variables:
    model -- the model containing the model_parts
    settings -- Kratos parameters containing solver settings.
    """

    def __init__(self, model, params):
        """ The default constructor of the class

        Keyword arguments:
        self -- It signifies an instance of a class.
        model -- the model containing the model_parts
        settings -- Kratos parameters containing solver settings.
        """
        KratosMultiphysics.Process.__init__(self)
        self.process = KratosMultiphysics.FromJSONCheckResultProcess(model, params)

    def ExecuteInitialize(self):
        """ This method is executed at the beginning to initialize the process

        Keyword arguments:
        self -- It signifies an instance of a class.
        """
        self.process.ExecuteInitialize()

    def ExecuteFinalizeSolutionStep(self):
        """ This method is executed in order to finalize the current step

        Here the dictionary containing the solution is filled

        Keyword arguments:
        self -- It signifies an instance of a class.
        """
        self.process.ExecuteFinalizeSolutionStep()
        self.assertTrue(self.process.IsCorrectResult(), "Results do not coincide with the JSON reference\n" + self.process.GetErrorMessage())
