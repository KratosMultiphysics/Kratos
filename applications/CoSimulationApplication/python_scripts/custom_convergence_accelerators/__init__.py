# Adding the current directory to the path such that the modules can be imported with __import__ in the factory
from __future__ import print_function, absolute_import, division  # makes these scripts backward compatible with python 2.6 and 2.7
import sys, os
sys.path.append(os.path.dirname(__file__))

available_convergence_accelerators = {
    "constant_relaxation"  : "convergence_accelerator_constant",
    "aitken"               : "convergence_accelerator_aitken",
    "iqnils"               : "convergence_accelerator_iqnils",
    "mvqn"                 : "convergence_accelerator_mvqn"
    }

def CreateConvergenceAccelerator(settings, solvers, level):
    """This function creates and returns the convergence accelerator used for CoSimulation
    New convergence accelerators have to be registered by adding them to "available_convergence_accelerators"
    """
    if (type(settings) != dict):
        raise Exception("Input is expected to be provided as a python dictionary")

    accelerator_type = settings["type"]

    if accelerator_type in available_convergence_accelerators:
        accelerator_module = __import__(available_convergence_accelerators[accelerator_type])
        return accelerator_module.Create(settings, solvers, level)
    else:
        err_msg  = 'The requested convergence accelerator "' + accelerator_type + '" is not available!\n'
        err_msg += 'The available convergence accelerators are:\n'
        for available_accelerator in available_convergence_accelerators:
            err_msg += "\t" + available_accelerator + "\n"
        raise NameError(err_msg)