# Adding the current directory to the path such that the modules can be imported with __import__ in the factory
from __future__ import print_function, absolute_import, division  # makes these scripts backward compatible with python 2.6 and 2.7
import sys, os
current_dir_name = os.path.dirname(__file__)
sys.path.append(current_dir_name)

# Append all the sub directories containing the solver interfaces to the path
for dirName, subdirList, fileList in os.walk(current_dir_name):
    sys.path.append(dirName)

available_solver_interfaces = {
    "kratos"  : "convergence_accelerator_constant",
    "dummy"               : "convergence_accelerator_aitken",
    "su2"               : "convergence_accelerator_iqnils",
    ""                 : "convergence_accelerator_mvqn"
    }

def CreateSolverInterface(settings, solvers, level):
    """
        This function creates and returns the solver interface used for CoSimulation
        New solver interface have to be registered by adding them to "available_convergence_accelerators"
    """
    if (type(settings) != dict):
        raise Exception("Input is expected to be provided as a python dictionary")

    solver_interface_type = settings["type"]

    if solver_interface_type in available_solver_interfaces:
        solver_interface_module = __import__(available_solver_interfaces[solver_interface_type])
        return solver_interface_module.Create(settings, solvers, level)
    else:
        err_msg  = 'The requested solver interface "' + solver_interface_type + '" is not available!\n'
        err_msg += 'Available solver interfaces are :\n'
        for avail_accelerator in available_convergence_accelerators:
            err_msg += "\t" + avail_accelerator + "\n"
        raise NameError(err_msg)