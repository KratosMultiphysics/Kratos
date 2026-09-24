# Importing the Kratos Library
import KratosMultiphysics as KM

from KratosMultiphysics import python_linear_solver_factory as linear_solver_factory
import KratosMultiphysics.kratos_utilities as kratos_utils

from typing import Union

def ConstructSolver(settings: KM.Parameters) -> Union[KM.LinearSolver, KM.ComplexLinearSolver]:
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

    else:
        return linear_solver_factory.ConstructSolver(settings)

