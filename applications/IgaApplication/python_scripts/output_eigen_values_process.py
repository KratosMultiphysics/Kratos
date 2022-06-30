# Importing the Kratos Library
import KratosMultiphysics
import KratosMultiphysics.IgaApplication as IGA

def Factory(settings, model):
    if not (isinstance(settings, KratosMultiphysics.Parameters)):
        raise Exception("Expected input shall be a Parameters object, encapsulating a json string")
    return IgaEigenValuesOutputProcess(model, settings["Parameters"])

class IgaEigenValuesOutputProcess(KratosMultiphysics.Process):

    def __init__(self, model, params):
        KratosMultiphysics.Process.__init__(self)
        print(params)
        self.process = IGA.OutputEigenValuesProcess(model, params)

    def ExecuteFinalize(self):
        print("Finalizeeee")
        self.process.ExecuteFinalize()