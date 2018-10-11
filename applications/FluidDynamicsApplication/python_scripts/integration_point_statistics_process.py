# Importing the Kratos Library
import KratosMultiphysics

# Check that applications were imported in the main script
KratosMultiphysics.CheckRegisteredApplications("FluidDynamicsApplication")

# Import applications
import KratosMultiphysics.FluidDynamicsApplication as KratosCFD

def Factory(settings, model):
    if(type(settings) != KratosMultiphysics.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")

    # TODO: pass the model directly
    model_part = model.GetModelPart(settings["model_part_name"].GetString())

    return KratosCFD.IntegrationPointStatisticsProcess(model_part, settings["Parameters"])