import KratosMultiphysics as Kratos
import KratosMultiphysics.FluidDynamicsApplication as KratosCFD

from KratosMultiphysics import IsDistributedRun
from KratosMultiphysics.kratos_utilities import CheckIfApplicationsAvailable

if (IsDistributedRun()
        and CheckIfApplicationsAvailable("TrilinosApplication")):
    import KratosMultiphysics.TrilinosApplication as KratosTrilinos


def TrilinosPeriodicBlockBuilderAndSolver(linear_solver, communicator):
    return KratosTrilinos.TrilinosBlockBuilderAndSolverPeriodic(
        communicator, 30, linear_solver,
        KratosCFD.PATCH_INDEX)


def TrilinosBlockBuilderAndSolver(linear_solver, communicator):
    return KratosTrilinos.TrilinosBlockBuilderAndSolver(
        communicator, 30, linear_solver)


def PeriodicBlockBuilderAndSolver(linear_solver, _):
    return KratosCFD.ResidualBasedBlockBuilderAndSolverPeriodic(
        linear_solver, KratosCFD.PATCH_INDEX)


def BlockBuilderAndSolver(linear_solver, _):
    return Kratos.ResidualBasedBlockBuilderAndSolver(linear_solver)