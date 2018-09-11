# Adding the current directory to the path such that the modules can be imported with __import__ in the factory
from __future__ import print_function, absolute_import, division  # makes these scripts backward compatible with python 2.6 and 2.7
import sys, os
current_dir_name = os.path.dirname(__file__)
sys.path.append(current_dir_name)

# Append all the sub directories containing the solver interfaces to the path
for dirName, subdirList, fileList in os.walk(current_dir_name):
    sys.path.append(dirName)

"""
This is a map of the name of the available solver interfaces to be specified in
JSON file and their python module (file) name. New solver interfaces should be registered
here with an additional entry.
eg : "name_in_JSON" : "python module(file) name"
"""
available_solver_interfaces = {
    "kratos_structure"      : "kratos_structure_co_simulation_solver_wrapper",
    "kratos_fluid"          : "kratos_fluid_co_simulation_solver_wrapper",
    "dummy"                 : "dummy_co_simulation_solver_wrapper",
    "su2"                   : "su2_co_simulation_solver_wrapper",
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





"""
This is a map of the name of the available IO types to be specified in
JSON file and their python module (file) name.
New IOs should be registered here with an additional entry.
eg : "name_in_JSON" : "python module(file) name"
"""
available_ios = {
    "kratos" : "kratos_io",
    "sdof"   : "su2_io",
}

def CreateIO(io_name, solvers, solver_name, level):
    """
    This function creates and returns the IO used for CoSimulation
    New IOs have to be registered by adding them to "available_ios"
    """
    if io_type in available_ios:
        io_module = __import__(available_ios[io_type])
        return io_module.Create(solvers, solver_name, level)
    else:
        err_msg  = 'The requested IO "' + io_name + '" is not available!\n'
        err_msg += 'The available IOs are :\n'
        for avail_io in available_ios:
            err_msg += "\t" + avail_io + "\n"
        raise NameError(err_msg)