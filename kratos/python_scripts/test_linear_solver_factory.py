from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
#
class bicg_solver_config:
    solver_type = "BiConjugate gradient stabilized"
    scaling = False
    # preconditioner_type = None
    max_iteration = 500
    tolerance = 1e-6

#


class cg_solver_config:
    solver_type = "Conjugate gradient"
    scaling = True
    preconditioner_type = "DiagonalPreconditioner"
    max_iteration = 1000
    tolerance = 1e-3

#


class gmres_solver_config:
    solver_type = "GMRES"
    scaling = True
    preconditioner_type = "None"
    max_iteration = 1000
    tolerance = 1e-3

#


class deflated_solver_config:
    solver_type = "Deflated Conjugate gradient"
    scaling = True
    max_iteration = 1000
    tolerance = 1e-3
    assume_constant_structure = True
    max_reduced_size = 1000

#


class amgcl_solver_config:
    solver_type = "AMGCL"
    scaling = False
    preconditioner_type = None
    max_iteration = 500
    tolerance = 1e-6
    smoother_type = "ILU0"
    krylov_type = "GMRES"

#


class superlu_config:  # note that i repeat the name of "velocity_linear_solver_config" which is not a problem since it is within a class
    solver_type = "Super LU"
    scaling = False

#


class superlu_iterative_solver_configure:
    solver_type = "SuperLUIterativeSolver"
    scaling = False
    tolerance = 1e-7
    max_iteration = 300
    gmres_krylov_space_dimension = 100
    DropTol = 1e-4
    FillTol = 1e-2
    ilu_level_of_fill = 5

#


class pastix_direct_solver_configure:
    solver_type = "PastixDirect"
    verbosity = 0
    is_symmetric = False

#


class pastix_iterative_solver_configure:
    solver_type = "PastixIterative"
    scaling = False
    tolerance = 1e-7
    max_iteration = 300
    gmres_krylov_space_dimension = 100
    ilu_level_of_fill = 5
    verbosity = 0
    is_symmetric = False


#
class block_UP_iterative_solver_configure:
    solver_type = "GMRES-UP Block"
    scaling = False
    tolerance = 1e-3
    max_iteration = 5
    gmres_krylov_space_dimension = 5

    class velocity_block_configuration:
        solver_type = "PastixIterative"
        scaling = False
        tolerance = 1e-6
        max_iteration = 300
        gmres_krylov_space_dimension = 100
        ilu_level_of_fill = 3
        verbosity = 0
        is_symmetric = False

    class pressure_block_configuration:
        solver_type = "AMGCL"
        scaling = False
        max_iteration = 50
        tolerance = 1e-3
        smoother_type = "DAMPED_JACOBI"
        krylov_type = "GMRES"

#
# KRATOS IS IMPORTED ONLY FROM THIS POINT ON
#
from KratosMultiphysics import *
from KratosMultiphysics.ExternalSolversApplication import *
import linear_solver_factory


linear_solver = linear_solver_factory.ConstructSolver(bicg_solver_config)

linear_solver = linear_solver_factory.ConstructSolver(cg_solver_config)

linear_solver = linear_solver_factory.ConstructSolver(gmres_solver_config)

linear_solver = linear_solver_factory.ConstructSolver(deflated_solver_config)

linear_solver = linear_solver_factory.ConstructSolver(amgcl_solver_config)

linear_solver = linear_solver_factory.ConstructSolver(
    superlu_iterative_solver_configure)

linear_solver = linear_solver_factory.ConstructSolver(
    pastix_direct_solver_configure)

linear_solver = linear_solver_factory.ConstructSolver(
    pastix_iterative_solver_configure)

linear_solver = linear_solver_factory.ConstructSolver(
    block_UP_iterative_solver_configure)
