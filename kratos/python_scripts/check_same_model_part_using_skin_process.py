# Importing the Kratos Library
import KratosMultiphysics

def Factory(settings, model):
    if not isinstance(settings, KratosMultiphysics.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return CheckSameModelPartUsingSkinDistanceProcess(model, settings["Parameters"])

class CheckSameModelPartUsingSkinDistanceProcess(KratosMultiphysics.Process):
    """This process checks that the model part is the same using the skin distance

    Only the member variables listed below should be accessed directly.

    Public member variables:
    Model -- the container of the different model parts.
    settings -- Kratos parameters containing solver settings.
    """

    def __init__(self, model, settings):
        """ The default constructor of the class

        Keyword arguments:
        self -- It signifies an instance of a class.
        model -- the container of the different model parts.
        settings -- Kratos parameters containing solver settings.
        """
        KratosMultiphysics.Logger.PrintWarning("DEPRECATION-Warning", "This python process will be removed in the future. Please use the C++ version instead with the register instead")
        KratosMultiphysics.Process.__init__(self)

        # Assigning values
        self.model = model

        # Get first model part (must be defined anyway)
        self.model_part = self.model[settings["skin_model_part_1_name"].GetString()]
        if not self.model_part.ProcessInfo.Has(KratosMultiphysics.DOMAIN_SIZE):
            raise ValueError("DOMAIN_SIZE nor defined in ProcessInfo")


        # Create the process
        if self.dimension == 2:
            self.process = KratosMultiphysics.CheckSameModelPartUsingSkinDistanceProcess2D(self.model, settings)
        else:
            self.process = KratosMultiphysics.CheckSameModelPartUsingSkinDistanceProcess3D(self.model, settings)

    def ExecuteBeforeSolutionLoop(self):
        """This method is executed in before initialize the solution step

        Keyword arguments:
        self -- It signifies an instance of a class.
        """
        self.process.Execute()
