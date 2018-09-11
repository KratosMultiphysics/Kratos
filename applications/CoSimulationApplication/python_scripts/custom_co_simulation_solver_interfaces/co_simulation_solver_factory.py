from __future__ import print_function, absolute_import, division  # makes these scripts backward compatible with python 2.6 and 2.7
"""
This is a map of the name of available convergence accelerators to be specified in
JSON file and their python module (file) name. New accelerators should be registered here by an
additional entry
eg : "name_in_JSON" : "python module(file) name"
"""
available_solver_interfaces = {
    "kratos_structural"    : "kratos_structural_co_simulation_solver",
    "kratos_fluid"         : "kratos_fluid_co_simulation_solver",
    "dummy"                : "dummy_co_simulation_solver"
    }

def CreateSolverInterface(settings):
    """
    This function creates and returns the convergence accelerator used for CoSimulation
    New convergence accelerators have to be registered by adding them to "available_convergence_accelerators"
    """
    if (type(settings) != dict):
        raise Exception("Input is expected to be provided as a python dictionary")

    solver_type = settings["solver_type"]

    if solver_type in available_solver_interfaces:
        solvers_module_name = "custom_co_simulation_solver_interfaces"
        __import__(solvers_module_name)
        solver_module = __import__(available_solver_interfaces[solver_type])
        return solver_module.Create(settings)
    else:
        err_msg  = 'The requested solver interface "' + solver_type + '" is not available!\n'
        err_msg += 'The available solver interfaces are:\n'
        for available_solver_interface in available_solver_interfaces:
            err_msg += "\t" + available_solver_interface + "\n"
        raise NameError(err_msg)