import sys
import importlib

import KratosMultiphysics
from KratosMultiphysics import kratos_utilities
from KratosMultiphysics.RomApplication import rom_solver


def _GetAvailableSolverWrapperModules():
    return {
        "KratosMultiphysics.FluidDynamicsApplication"       : "python_solvers_wrapper_fluid",
        "KratosMultiphysics.StructuralMechanicsApplication" : "python_solvers_wrapper_structural",
        "KratosMultiphysics.ConvectionDiffusionApplication" : "python_solvers_wrapper_convection_diffusion",
        "KratosMultiphysics.CompressiblePotentialFlowApplication" : "python_solvers_wrapper_compressible_potential"
    }


def CreateSolverByParameters(model, solver_settings, parallelism, analysis_stage_module_name):

    if not isinstance(model, KratosMultiphysics.Model):
        raise Exception("input is expected to be provided as a Kratos Model object")

    if not isinstance(solver_settings, KratosMultiphysics.Parameters):
        raise Exception("input is expected to be provided as a Kratos Parameters object")

    # Get the corresponding application from the analysis_stage path
    split_analysis_stage_module_name = analysis_stage_module_name.split('.')
    application_module_name = split_analysis_stage_module_name[0] + "." + split_analysis_stage_module_name[1]
    if not kratos_utilities.CheckIfApplicationsAvailable(split_analysis_stage_module_name[1]):
        raise Exception("Module {} is not available.".format(application_module_name))

    # Filter and retrieve the Python solvers wrapper from the corresponding application
    #TODO: This filtering wouldn't be required if we were using a unified solvers wrapper module name
    available_modules = _GetAvailableSolverWrapperModules()

    if application_module_name in available_modules:
        solvers_wrapper_module_module_name = available_modules[application_module_name]
    else:
        err_msg = "Python module \'{0}\' is not available. Make sure \'{1}\' is compiled and implemented.\n".format(
            application_module_name, split_analysis_stage_module_name[1])
        err_msg += "Currently implemented applications are:\n"
        err_msg += "".join(" - {}\n".format(key) for key in available_modules)
        err_msg += "To add a new implementation, do so in '{}' in {}".format(
            _GetAvailableSolverWrapperModules.__name__, __file__)
        raise Exception(err_msg)
    solvers_wrapper_module = importlib.import_module(application_module_name + "." + solvers_wrapper_module_module_name)

    # Create a prototype class instance and get the module and name of the solver to be used as base
    # Note that an auxiliary Kratos parameter settings without the rom_settings field is created to avoid the defaults error thrown
    # Note that an auxiliary Kratos model is also created to avoid creating the main_model_part in the prototype class instance
    #TODO: We could do the same exercise as we do in the stage (module_name to ClassName equal to ModuleName if we standarize the solver names)
    aux_solver_settings = solver_settings.Clone()
    aux_solver_settings.RemoveValue("rom_settings")
    aux_solver_settings.RemoveValue("projection_strategy")
    aux_solver_settings.RemoveValue("assembling_strategy")
    aux_base_solver_instance = solvers_wrapper_module.CreateSolverByParameters(KratosMultiphysics.Model(), aux_solver_settings, parallelism)

    # Create the ROM solver from the base solver
    rom_solver_instance = rom_solver.CreateSolver(type(aux_base_solver_instance), model, solver_settings)

    return rom_solver_instance

def CreateSolver(model, custom_settings):

    if (type(model) != KratosMultiphysics.Model):
        raise Exception("input is expected to be provided as a Kratos Model object")

    if (type(custom_settings) != KratosMultiphysics.Parameters):
        raise Exception("input is expected to be provided as a Kratos Parameters object")

    parallelism = custom_settings["problem_data"]["parallel_type"].GetString()
    analysis_stage = custom_settings["analysis_stage"].GetString()
    solver_settings = custom_settings["solver_settings"]

    return CreateSolverByParameters(model, solver_settings, parallelism, analysis_stage)

