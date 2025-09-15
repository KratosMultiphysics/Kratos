# Importing the Kratos Library
import KratosMultiphysics
import KratosMultiphysics.IgaApplication as IGA
from KratosMultiphysics import kratos_utilities

def Factory(settings, model):
    if not isinstance(settings, KratosMultiphysics.Parameters):
        raise Exception("Expected input shall be a Parameters object, encapsulating a json string")
    return MapNurbsVolumeResultsToEmbeddedGeometryProcess(model, settings["Parameters"])

class MapNurbsVolumeResultsToEmbeddedGeometryProcess(KratosMultiphysics.Process):
    def __init__(self, model, params):
        KratosMultiphysics.Process.__init__(self)
        self.params = params
        self.process = IGA.MapNurbsVolumeResultsToEmbeddedGeometryProcess(model, params)

    def ExecuteBeforeOutputStep(self):
        self.process.MapVariables()
