from KratosMultiphysics.CoSimulationApplication.factories import base_factory

def CreateIO(io_settings, model, solver_name, data_comm, io_name):
    """This function creates and returns the IO used for CoSimulation"""
    return base_factory.Create(io_settings, [model, solver_name, data_comm], "KratosMultiphysics.CoSimulationApplication.solver_wrappers", io_name)
