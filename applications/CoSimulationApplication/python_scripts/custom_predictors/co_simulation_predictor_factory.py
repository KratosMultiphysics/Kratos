from __future__ import print_function, absolute_import, division  # makes these scripts backward compatible with python 2.6 and 2.7
"""
This is a map of the name of the available predictor types to be specified in
JSON file and their python module (file) name.
New IOs should be registered here with an additional entry.
eg : "name_in_JSON" : "python module(file) name"
"""
available_predictors = {
    "linear"  : "linear_predictor",
}

def CreatePredictor(settings, solver):
    """
    This function creates and returns the predictor used for CoSimulation
    New predictors have to be registered by adding them to "available_predictors" above
    """
    predictor_type = settings["type"].GetString()
    if predictor_type in available_predictors:
        module_name = available_predictors[predictor_type]
        module_full = "custom_predictors."+module_name
        accelerator_module = __import__(module_full,fromlist=[module_name])
        return accelerator_module.Create(settings, solver)
    else:
        err_msg  = 'The requested predictor "' + predictor_type + '" is not available!\n'
        err_msg += 'The available predictors are:\n'
        for avail_predictor in available_predictors:
            err_msg += "\t" + avail_predictor + "\n"
        raise NameError(err_msg)