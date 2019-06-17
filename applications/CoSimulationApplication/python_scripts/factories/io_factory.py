from __future__ import print_function, absolute_import, division  # makes these scripts backward compatible with python 2.6 and 2.7

def CreateIO(io_name, model, io_settings):
    """This function creates and returns the IO used for CoSimulation
    The io-module has to be on the PYTHONPATH
    Naming-Convention: The module-file has to end with "_io"
    """

    module_name = io_name

    # TODO come up with sth better, this is hardcoded to Kratos!
    module_full = "KratosMultiphysics.CoSimulationApplication.solver_wrappers."+module_name

    io_module = __import__(module_full,fromlist=[module_name])
    return io_module.Create(model, io_settings)
