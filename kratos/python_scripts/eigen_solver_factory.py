from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7

# Importing the Kratos Library
import KratosMultiphysics as KM

from KratosMultiphysics import python_linear_solver_factory as linear_solver_factory
import KratosMultiphysics.kratos_utilities as kratos_utils

if kratos_utils.CheckIfApplicationsAvailable("EigenSolversApplication"):
    import KratosMultiphysics.EigenSolversApplication as EiSA

if kratos_utils.CheckIfApplicationsAvailable("ExternalSolversApplication"):
    import KratosMultiphysics.ExternalSolversApplication as ExSA

def ConstructSolver(settings):
    if not isinstance(settings, KM.Parameters):
        raise Exception("Input is expected to be provided as a Kratos Parameters object")

    solver_type = settings["solver_type"].GetString()

    if solver_type == "eigen_eigensystem":
        if kratos_utils.CheckIfApplicationsAvailable("EigenSolversApplication"):
            eigen_solver = EiSA.EigensystemSolver(settings)
            return eigen_solver
        else:
            raise Exception("EigenSolversApplication not available")

    linear_solver = CreateLinearSolver(settings)

    if solver_type == "power_iteration_eigenvalue_solver":
        eigen_solver = KM.PowerIterationEigenvalueSolver( settings, linear_solver)
    elif solver_type == "power_iteration_highest_eigenvalue_solver":
        eigen_solver = KM.PowerIterationHighestEigenvalueSolver( settings, linear_solver)
    elif solver_type == "rayleigh_quotient_iteration_eigenvalue_solver":
        eigen_solver = KM.RayleighQuotientIterationEigenvalueSolver( settings, linear_solver)
    elif solver_type == "FEAST" or solver_type == "feast":
        if kratos_utils.CheckIfApplicationsAvailable("ExternalSolversApplication"):
            eigen_solver = ExSA.FEASTSolver(settings, linear_solver)
        else:
            raise Exception("ExternalSolversApplication not available")
    else:
        raise Exception("Solver type not found. Asking for :" + solver_type)

    return eigen_solver

def CreateLinearSolver(settings):
    if not settings.Has("linear_solver_settings"):
        settings.AddEmptyValue("linear_solver_settings")
    linear_solver_configuration = settings["linear_solver_settings"]
    if linear_solver_configuration.Has("solver_type"): # user specified a linear solver
        return linear_solver_factory.ConstructSolver(linear_solver_configuration)
    else:
        # Using a default linear solver (selecting the fastest one available)
        linear_solvers_by_speed = [
            "pardiso_lu", # EigenSolversApplication (if compiled with Intel-support)
            "sparse_lu",  # EigenSolversApplication
            "pastix",     # ExternalSolversApplication (if Pastix is included in compilation)
            "super_lu",   # ExternalSolversApplication
            "skyline_lu_factorization" # in Core, always available, but slow
        ]

        for solver_name in linear_solvers_by_speed:
            if KM.LinearSolverFactory().Has(solver_name):
                linear_solver_configuration.AddEmptyValue("solver_type").SetString(solver_name)
                KM.Logger.PrintInfo('::[MechanicalSolver]:: ',\
                    'Using "' + solver_name + '" as default linear solver')
                return KM.LinearSolverFactory().Create(linear_solver_configuration)

    raise Exception("Linear-Solver could not be constructed!")
