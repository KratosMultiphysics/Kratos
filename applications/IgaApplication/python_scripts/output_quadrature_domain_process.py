# Importing the Kratos Library
import KratosMultiphysics
import KratosMultiphysics.IgaApplication as IGA


def Factory(settings, model):
    if(type(settings) != KratosMultiphysics.Parameters):
        raise Exception("Expected input shall be a Parameters object, encapsulating a json string")
    return IgaOutputIntegrationDomainProcess(model, settings["Parameters"])

class IgaOutputIntegrationDomainProcess(KratosMultiphysics.Process):

    def __init__(self, model, params):
        KratosMultiphysics.Process.__init__(self)
        self.process = IGA.OutputQuadratureDomainProcess(model, params)

    def ExecuteBeforeSolutionLoop(self):
        self.process.ExecuteBeforeSolutionLoop()