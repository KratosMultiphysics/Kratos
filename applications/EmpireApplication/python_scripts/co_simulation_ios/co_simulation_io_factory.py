from __future__ import print_function, absolute_import, division  # makes these scripts backward compatible with python 2.6 and 2.7

def CreateIO(io_name, solvers, solver_name, level):
    """This function creates and returns the IO used for CoSimulation
    The io-module has to be on the PYTHONPATH
    Naming-Convention: The module-file has to end with "_io"
    """

    io_module_name = io_name + "_io"

    io_module = __import__(io_module_name)
    return io_module.Create(solvers, solver_name, level)
