# Importing the Kratos Library
import KratosMultiphysics as KM

from KratosMultiphysics import python_linear_solver_factory as linear_solver_factory
import KratosMultiphysics.kratos_utilities as kratos_utils

def ConstructSolver(settings):
    if not isinstance(settings, KM.Parameters):
        raise Exception("Input is expected to be provided as a Kratos Parameters object")

    solver_type = settings["solver_type"].GetString()

    if solver_type == "eigen_eigensystem":
        if kratos_utils.CheckIfApplicationsAvailable("LinearSolversApplication"):
            from KratosMultiphysics import LinearSolversApplication
            eigen_solver = LinearSolversApplication.EigensystemSolver(settings)
            return eigen_solver
        else:
            raise Exception("LinearSolversApplication not available")

    if solver_type == "spectra_sym_g_eigs_shift":
        if kratos_utils.CheckIfApplicationsAvailable("LinearSolversApplication"):
            from KratosMultiphysics import LinearSolversApplication
            eigen_solver = LinearSolversApplication.SpectraSymGEigsShiftSolver(settings)
            return eigen_solver
        else:
            raise Exception("LinearSolversApplication not available")
    if solver_type == "spectra_g_eigs_shift":
        if kratos_utils.CheckIfApplicationsAvailable("LinearSolversApplication"):
            from KratosMultiphysics import LinearSolversApplication
            eigen_solver = LinearSolversApplication.SpectraGEigsShiftSolver(settings)
            return eigen_solver
        else:
            raise Exception("LinearSolversApplication not available")


    elif solver_type == "dense_eigensolver":
        if kratos_utils.CheckIfApplicationsAvailable("LinearSolversApplication"):
            from KratosMultiphysics import LinearSolversApplication
            eigen_solver = LinearSolversApplication.DenseEigenvalueSolver(settings)
            return eigen_solver
        else:
            raise Exception("LinearSolversApplication not available")

    elif solver_type == "feast":
        if kratos_utils.CheckIfApplicationsAvailable("LinearSolversApplication"):
            from KratosMultiphysics import LinearSolversApplication
            if LinearSolversApplication.HasFEAST():
                is_symmetric = settings["symmetric"].GetBool() if settings.Has("symmetric") else True
                eigen_solver = LinearSolversApplication.FEASTSymmetricEigensystemSolver(settings) if is_symmetric else LinearSolversApplication.FEASTGeneralEigensystemSolver(settings)
                return eigen_solver
            else:
                raise Exception("FEAST not available in LinearSolversApplication")
        else:
            raise Exception("LinearSolversApplication not available")

    elif solver_type == "feast_complex":
        if kratos_utils.CheckIfApplicationsAvailable("LinearSolversApplication"):
            from KratosMultiphysics import LinearSolversApplication
            if LinearSolversApplication.HasFEAST():
                is_symmetric = settings["symmetric"].GetBool() if settings.Has("symmetric") else True
                eigen_solver = LinearSolversApplication.ComplexFEASTSymmetricEigensystemSolver(settings) if is_symmetric else LinearSolversApplication.ComplexFEASTGeneralEigensystemSolver(settings)
                return eigen_solver
            else:
                raise Exception("FEAST not available in LinearSolversApplication")
        else:
            raise Exception("LinearSolversApplication not available")

    linear_solver_configuration = settings["linear_solver_settings"]
    if linear_solver_configuration.Has("solver_type"): # user specified a linear solver
        linear_solver = linear_solver_factory.ConstructSolver(linear_solver_configuration)
    else:
        linear_solver = linear_solver_factory.CreateFastestAvailableDirectLinearSolver()

    if solver_type == "power_iteration_eigenvalue_solver":
        eigen_solver = KM.PowerIterationEigenvalueSolver( settings, linear_solver)
    elif solver_type == "power_iteration_highest_eigenvalue_solver":
        eigen_solver = KM.PowerIterationHighestEigenvalueSolver( settings, linear_solver)
    elif solver_type == "rayleigh_quotient_iteration_eigenvalue_solver":
        eigen_solver = KM.RayleighQuotientIterationEigenvalueSolver( settings, linear_solver)
    else:
        raise Exception("Solver type not found. Asking for :" + solver_type)

    return eigen_solver
