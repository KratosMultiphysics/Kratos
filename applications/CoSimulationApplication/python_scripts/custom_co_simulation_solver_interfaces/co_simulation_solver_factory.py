from __future__ import print_function, absolute_import, division  # makes these scripts backward compatible with python 2.6 and 2.7
"""
This is a map of the name of available convergence accelerators to be specified in
JSON file and their python module (file) name. New accelerators should be registered here by an
additional entry
eg : "name_in_JSON" : "python module(file) name"
"""
available_solver_interfaces = {
    "kratos_structural"    : "kratos_interfaces.kratos_structural_co_simulation_solver",
    "kratos_fluid"         : "kratos_interfaces.kratos_fluid_co_simulation_solver",
    "dummy"                : "dummy_co_simulation_solver"
    }

def CreateSolverInterface(solver_name, settings):
    """
    This function creates and returns the convergence accelerator used for CoSimulation
    New convergence accelerators have to be registered by adding them to "available_convergence_accelerators"
    """

    solver_type = settings["solver_type"].GetString()

    if solver_type in available_solver_interfaces:
        solver_module_name = available_solver_interfaces[solver_type]
        module_full = 'custom_co_simulation_solver_interfaces.'+solver_module_name
        solver_module = __import__(module_full,fromlist=[solver_module_name])
        return solver_module.Create(solver_name, settings)
    else:
        err_msg  = 'The requested solver interface "' + solver_type + '" is not available!\n'
        err_msg += 'The available solver interfaces are:\n'
        for available_solver_interface in available_solver_interfaces:
            err_msg += "\t" + available_solver_interface + "\n"
        raise NameError(err_msg)
