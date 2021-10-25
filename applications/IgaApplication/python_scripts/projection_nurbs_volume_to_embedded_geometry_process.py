# Importing the Kratos Library
import KratosMultiphysics
import KratosMultiphysics.IgaApplication as IGA


def Factory(settings, model):
    if not isinstance(settings, KratosMultiphysics.Parameters):
        raise Exception("Expected input shall be a Parameters object, encapsulating a json string")
    return ProjectionNurbsVolumeToEmbeddedGeometryProcess(model, settings["Parameters"])

class ProjectionNurbsVolumeToEmbeddedGeometryProcess(KratosMultiphysics.Process):

    def __init__(self, model, params):
        KratosMultiphysics.Process.__init__(self)
        self.process = IGA.ProjectionNurbsVolumeToEmbeddedGeometryProcess(model, params)

    def ExecuteBeforeOutputStep(self):
        self.process.ExecuteBeforeOutputStep()