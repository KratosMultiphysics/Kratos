# Importing the Kratos Library
import KratosMultiphysics
import KratosMultiphysics.IgaApplication as IGA


def Factory(settings, model):
    if not isinstance(settings, KratosMultiphysics.Parameters):
        raise Exception("Expected input shall be a Parameters object, encapsulating a json string")
    return OutputIntegrationDomainProcess(model, settings["Parameters"])

class OutputIntegrationDomainProcess(KratosMultiphysics.Process):

    def __init__(self, model, params):
        KratosMultiphysics.Process.__init__(self)
        self.process = IGA.OutputQuadratureDomainProcess(model, params)

    def ExecuteBeforeSolutionLoop(self):
        self.process.ExecuteBeforeSolutionLoop()