import KratosMultiphysics
import KratosMultiphysics.LinearSolversApplication as KratosLSA

def Factory(settings, model):
    if (not isinstance(model, KratosMultiphysics.Model)):
        raise Exception(
            "expected input shall be a Model object, encapsulating a json string"
        )
    if (not isinstance(settings, KratosMultiphysics.Parameters)):
        raise Exception(
            "expected input shall be a Parameters object, encapsulating a json string"
        )
    return KratosLSA.MaxSingularValueDecompositionProcess(model, settings["Parameters"])

