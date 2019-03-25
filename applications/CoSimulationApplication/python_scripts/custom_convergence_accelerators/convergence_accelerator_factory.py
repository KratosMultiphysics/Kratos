from __future__ import print_function, absolute_import, division  # makes these scripts backward compatible with python 2.6 and 2.7


def CreateConvergenceAccelerator(settings, solver):
    """
    This function creates and returns the convergence accelerator.
    """
    accelerator_type = settings["type"]
    module_name = accelerator_type.GetString()
    module_full = "KratosMultiphysics.CoSimulationApplication.custom_convergence_accelerators.convergence_accelerator_"+module_name
    accelerator_module = __import__(module_full,fromlist=[module_name])
    return accelerator_module.Create(settings, solver)
