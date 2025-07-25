# Importing the Kratos Library
import KratosMultiphysics
import KratosMultiphysics.IgaApplication as IGA


def Factory(settings, model):
    if not isinstance(settings, KratosMultiphysics.Parameters):
        raise Exception("Expected input shall be a Parameters object, encapsulating a json string")
    return ComputeBeamVectorsProcess(model, settings["Parameters"])

class ComputeBeamVectorsProcess(KratosMultiphysics.Process):
    """This process computes the tangential (T0) and normal (N0) vectors for isogeometric beam elements during initialization.
    """
    def __init__(self, model, params):
        KratosMultiphysics.Process.__init__(self)
        self.process = IGA.ComputeBeamVectorsProcess(model, params)

    def ExecuteInitialize(self):
        self.process.ExecuteInitialize()     