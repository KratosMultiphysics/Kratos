from __future__ import print_function, absolute_import, division  # makes these scripts backward compatible with python 2.6 and 2.7
"""
This is a map of the name of the available extrapolator types to be specified in
JSON file and their python module (file) name.
New IOs should be registered here with an additional entry.
eg : "name_in_JSON" : "python module(file) name"
"""
available_extrapolators = {
    "linear"  : "linear_extrapolator",
}

def CreateExtrapolators(settings, solver):
    """
    This function creates and returns the extrapolator used for CoSimulation
    New extrapolators have to be registered by adding them to "available_extrapolators" above
    """
    accelerator_type = settings["type"].GetString()
    if accelerator_type in available_convergence_accelerators:
        module_name = available_convergence_accelerators[accelerator_type]
        module_full = "custom_extrapolators."+module_name
        accelerator_module = __import__(module_full,fromlist=[module_name])
        return accelerator_module.Create(settings, solver)
    else:
        err_msg  = 'The requested extrapolator "' + accelerator_type + '" is not available!\n'
        err_msg += 'The available extrapolators are:\n'
        for avail_accelerator in available_convergence_accelerators:
            err_msg += "\t" + avail_accelerator + "\n"
        raise NameError(err_msg)