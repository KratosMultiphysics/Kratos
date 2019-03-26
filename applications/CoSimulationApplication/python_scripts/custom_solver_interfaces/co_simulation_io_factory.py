from __future__ import print_function, absolute_import, division  # makes backward compatible with python 2.6 and 2.7


def CreateIO(model, settings):
    """
    This function creates and returns the IO used for CoSimulation.
    """
    io_type = settings["io_type"]
    module_name = io_type.GetString()
    module_full = 'KratosMultiphysics.CoSimulationApplication.custom_solver_interfaces.'+module_name
    io_module = __import__(module_full, fromlist=[module_name])
    return io_module.Create(model, settings)
