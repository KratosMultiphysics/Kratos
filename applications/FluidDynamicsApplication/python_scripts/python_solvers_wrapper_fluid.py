from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7

import KratosMultiphysics

def CreateSolverByParameters(main_model_part, solver_settings, parallelism):
    
    solver_type = solver_settings["solver_type"].GetString()
    
    # Solvers for OpenMP parallelism
    if (parallelism == "OpenMP"):
        if (solver_type == "Monolithic"):
            solver_module_name = "navier_stokes_solver_vmsmonolithic"

        elif (solver_type == "FractionalStep"):
            solver_module_name = "navier_stokes_solver_fractionalstep"

        elif ((solver_type == "Embedded") or (solver_type == "EmbeddedDevelopment")):
            solver_module_name = "navier_stokes_embedded_solver"

        elif (solver_type == "EmbeddedAusas"):
            solver_module_name = "navier_stokes_embedded_ausas_solver"

        elif (solver_type == "Compressible"):
            solver_module_name = "navier_stokes_compressible_solver"

        elif (solver_type == "ThermallyCoupled"):
            solver_module_name = "coupled_fluid_thermal_solver"

        else:
            raise Exception("the requested solver type is not in the python solvers wrapper. Solver type is : "+solver_type)

    # Solvers for MPI parallelism
    elif (parallelism == "MPI"):
        if (solver_type == "Monolithic"):
            solver_module_name = "trilinos_navier_stokes_solver_vmsmonolithic"

        elif (solver_type == "FractionalStep"):
            solver_module_name = "trilinos_navier_stokes_solver_fractionalstep"

        elif (solver_type == "Embedded"):
            solver_module_name = "trilinos_navier_stokes_embedded_solver"

        elif (solver_type == "EmbeddedAusas"):
            solver_module_name = "trilinos_navier_stokes_embedded_ausas_solver"
  
        else:
            raise Exception("the requested solver type is not in the python solvers wrapper")
    else:
        raise Exception("parallelism is neither OpenMP nor MPI")

    solver_module = __import__(solver_module_name)
    solver = solver_module.CreateSolver(main_model_part, solver_settings)

    return solver

def CreateSolver(main_model_part, custom_settings):
    if (type(main_model_part) != KratosMultiphysics.ModelPart):
        raise Exception("input is expected to be provided as a Kratos ModelPart object")

    if (type(custom_settings) != KratosMultiphysics.Parameters):
        raise Exception("input is expected to be provided as a Kratos Parameters object")

    parallelism = custom_settings["problem_data"]["parallel_type"].GetString()
    solver_settings = custom_settings["solver_settings"]   
    
    return CreateSolverByParameters(main_model_part, solver_settings, parallelism)
