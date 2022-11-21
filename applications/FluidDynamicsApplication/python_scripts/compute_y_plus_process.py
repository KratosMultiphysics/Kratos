# Importing the Kratos Library
import KratosMultiphysics as Kratos
import KratosMultiphysics.FluidDynamicsApplication as KratosCFD

def Factory(settings, model):
    if not isinstance(settings, Kratos.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    if not isinstance(model, Kratos.Model):
        raise Exception("expected input shall be a model object")

    return KratosCFD.ComputeYPlusProcess(model, settings["Parameters"])




