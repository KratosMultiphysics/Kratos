# Importing the Kratos Library
import KratosMultiphysics as KM

# CoSimulation imports
import KratosMultiphysics.CoSimulationApplication.factories.solver_wrapper_factory as solver_wrapper_factory

def Create(settings, models, solver_name):
    input_file_name = settings["input_file"].GetString()
    settings.RemoveValue("input_file")
    if not input_file_name.endswith(".json"):
        input_file_name += ".json"

    with open(input_file_name,'r') as parameter_file:
        existing_parameters = KM.Parameters(parameter_file.read())
    for key, val in existing_parameters["solver_settings"].items():
        settings.AddValue(key, val)

    return solver_wrapper_factory.CreateSolverWrapper(settings, models, solver_name)
