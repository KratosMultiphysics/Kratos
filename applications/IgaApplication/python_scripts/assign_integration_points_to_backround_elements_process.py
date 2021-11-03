# Importing the Kratos Library
import KratosMultiphysics
import KratosMultiphysics.IgaApplication as IGA


def Factory(settings, model):
    if not isinstance(settings, KratosMultiphysics.Parameters):
        raise Exception("Expected input shall be a Parameters object, encapsulating a json string")
    return AssignIntegrationPointsToBackgroundElementsProcess(model, settings["Parameters"])

class AssignIntegrationPointsToBackgroundElementsProcess(KratosMultiphysics.Process):
    """This class assigns the integration points of an embedded geometry to the elements of the background mesh.

    Only the member variables listed below should be accessed directly.

    Public member variables:
    model -- the container of the different model parts.
    params -- Kratos parameters containing the settings.
    """
    def __init__(self, model, params):
        """ The default constructor of the class

        Keyword arguments:
        self -- It signifies an instance of a class.
        model -- the container of the different model parts.
        params -- Kratos parameters containing solver settings.
        """
        KratosMultiphysics.Process.__init__(self)
        self.process = IGA.AssignIntegrationPointsToBackgroundElementsProcess(model, params)

    def ExecuteBeforeOutputStep(self):
        """AssignIntegrationPoints must not be called in Initialize(), as output_quadrature_domain must be called a priori!

        However, the .cpp-Implementation takes care that this process is only called once!
        """
        self.process.AssignIntegrationPoints()
