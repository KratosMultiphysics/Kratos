# Importing the Kratos Library
import KratosMultiphysics

# Import applications
import KratosMultiphysics.FluidDynamicsApplication as KratosCFD

def Factory(settings, model):
    if(type(settings) != KratosMultiphysics.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")

    if not settings["Parameters"].Has("model_part_name"):
        msg = "Provided settings argument does not contain the model part name.\n"
        msg += "Please provide it as the Parameters/model_part_name (string) argument."
        raise Exception(msg)

    return KratosCFD.IntegrationPointStatisticsProcess(model, settings["Parameters"])