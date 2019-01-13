# Adding the current directory to the path such that the modules can be imported with __import__ in the factory
import sys, os
sys.path.append(os.path.dirname(__file__))

#Comment this is not up to date

"""
This is a map of the name of available convergence accelerators to be specified in
JSON file and their python module (file) name. New accelerators should be registered here by an
additional entry
eg : "name_in_JSON" : "python module(file) name"
"""
available_convergence_criteria = {
    "relative"  : "convergence_criteria_relative",
    "absolute"  : "convergence_criteria_absolute",
    }

def CreateConvergenceCriteria(settings, solvers, level):
    """
    This function creates and returns the convergence criteria used for CoSimulation
    New convergence criteria have to be registered by adding them to "availabcle_convergence_criteria"
    """
    if (type(settings) != dict):
        raise Exception("Input is expected to be provided as a python dictionary")

    criteria_type = settings["type"]

    if criteria_type in available_convergence_criteria:
        criteria_module = __import__(available_convergence_accelerators[criteria_type])
        return criteria_module.Create(settings, solvers, level)
    else:
        err_msg  = 'The requested convergence criteria "' + criteria_type + '" is not available!\n'
        err_msg += 'The available convergence criteria are:\n'
        for available_criteria in available_convergence_criteria:
            err_msg += "\t" + available_criteria + "\n"
        raise NameError(err_msg)