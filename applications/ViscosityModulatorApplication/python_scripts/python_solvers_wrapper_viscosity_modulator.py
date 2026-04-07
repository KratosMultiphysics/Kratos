
import KratosMultiphysics
from importlib import import_module

def CreateSolverByParameters(model, solver_settings, parallelism):

    if (type(model) != KratosMultiphysics.Model):
        raise Exception("input is expected to be provided as a Kratos Model object")

    if (type(solver_settings) != KratosMultiphysics.Parameters):
        raise Exception("input is expected to be provided as a Kratos Parameters object")

    solver_type = solver_settings["solver_type"].GetString()

    # Solvers for OpenMP parallelism
    if (parallelism == "OpenMP"):
        # Transient solvers
        if (solver_type == "transient" or solver_type == "Transient"):
            # If not provided, set implicit time integration as default
            if not solver_settings.Has("time_integration_method"):
                KratosMultiphysics.Logger.PrintWarning("Time integration method was not provided. Setting \'implicit\' as default.")
                solver_settings.AddEmptyValue("time_integration_method").SetString("implicit")
            time_integration_method = solver_settings["time_integration_method"].GetString()
            # Check transient integration method
            if time_integration_method == "implicit":
                solver_module_name = "viscosity_modulator_transient_solver"
            else:
                err_msg =  "The requested time integration method {} is not in the Python solvers wrapper\n".format(time_integration_method)
                err_msg += "Available options are: \"implicit\""
                raise Exception(err_msg)
        # Steady solver
        elif (solver_type == "stationary" or solver_type == "Stationary"):
            solver_module_name = "viscosity_modulator_stationary_solver"
        # Coupled CFD-scalar solvers (volume coupling by Boussinesq approximation)
        elif (solver_type == "scalar_coupled" or solver_type == "ScalarCoupled"):
            solver_module_name = "coupled_fluid_scalar_solver"
        # Coupled CFD-scalar solver with fully implicit Boussinesq body-force coupling
        elif (solver_type == "scalar_coupled_boussinesq" or solver_type == "ScalarCoupledBoussinesq"):
            solver_module_name = "coupled_fluid_scalar_solver_boussinesq"
        # Wrong solver check
        else:
            err_msg =  "The requested solver type {} is not in the python solvers wrapper\n".format(solver_type)
            err_msg += "Available options are: \"transient\", \"stationary\", \"scalar_coupled\", \"scalar_coupled_boussinesq\""
            raise Exception(err_msg)

    # Solvers for MPI parallelism
    elif (parallelism == "MPI"):
        # Transient solvers
        if (solver_type == "transient" or solver_type == "Transient"):
            # If not provided, set implicit time integration as default
            if not solver_settings.Has("time_integration_method"):
                KratosMultiphysics.Logger.PrintWarning("Time integration method was not provided. Setting \'implicit\' as default.")
                solver_settings.AddEmptyValue("time_integration_method").SetString("implicit")
            time_integration_method = solver_settings["time_integration_method"].GetString()
            # Check transient integration method
            if time_integration_method == "implicit":
                solver_module_name = "viscosity_modulator_transient_solver"
            else:
                err_msg =  "The requested time integration method {} is not MPI available yet\n".format(time_integration_method)
                err_msg += "Available option is \"implicit\""
                raise Exception(err_msg)
        # Steady solver
        elif (solver_type == "stationary" or solver_type == "Stationary"):
            solver_module_name = "viscosity_modulator_stationary_solver"
        # Wrong solver check
        else:
            err_msg =  "The requested solver type {} is not MPI available yet\n".format(solver_type)
            err_msg += "Available options are: \"transient\", \"stationary\""
            raise Exception(err_msg)
    else:
        err_msg =  "The requested parallel type \"" + parallelism + "\" is not available!\n"
        err_msg += "Available options are: \"OpenMP\", \"MPI\""
        raise Exception(err_msg)

    module_full = 'KratosMultiphysics.ViscosityModulatorApplication.' + solver_module_name
    solver = import_module(module_full).CreateSolver(model, solver_settings)

    return solver

def CreateSolver(model, custom_settings):

    if (type(model) != KratosMultiphysics.Model):
        raise Exception("input is expected to be provided as a Kratos Model object")

    if (type(custom_settings) != KratosMultiphysics.Parameters):
        raise Exception("input is expected to be provided as a Kratos Parameters object")

    parallelism = custom_settings["problem_data"]["parallel_type"].GetString()
    solver_settings = custom_settings["solver_settings"]

    return CreateSolverByParameters(model, solver_settings, parallelism)
