from __future__ import print_function, absolute_import, division  # makes these scripts backward compatible with python 2.6 and 2.7

def CreateConvergenceAccelerator(convergence_accelerator_settings, solver):
    """This function creates and returns the Convergence Accelerator used for CoSimulation
    The convergence-accelerator-module has to be on the PYTHONPATH
    Naming-Convention: The module-file has to end with "_convergence_accelerator"
    """
    module_name = convergence_accelerator_settings["type"].GetString()
    module_full = "KratosMultiphysics.CoSimulationApplication.convergence_accelerators."+module_name
    module_full += "_convergence_accelerator"

    accelerator_module = __import__(module_full,fromlist=[module_name])
    return accelerator_module.Create(convergence_accelerator_settings, solver)
