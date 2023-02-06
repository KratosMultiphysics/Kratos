import KratosMultiphysics as Kratos
import KratosMultiphysics.OptimizationApplication as KratosOA

def Factory(settings, model):
    if (not isinstance(model, Kratos.Model)):
        raise Exception(
            "expected input shall be a Model object, encapsulating a json string"
        )
    if (not isinstance(settings, Kratos.Parameters)):
        raise Exception(
            "expected input shall be a Parameters object, encapsulating a json string"
        )
    return KratosOA.ExtractBoundaryNodesProcess(model, settings["Parameters"])
