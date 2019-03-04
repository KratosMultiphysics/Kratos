from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
from KratosMultiphysics import *
from KratosMultiphysics.TrilinosApplication import *


def MultilevelLinearSolver(tolerance, max_iterations):
    # settings for the iterative solver
    aztec_parameters = ParameterList()
    # conjugate gradient solver for symmetric positive-definite
    # systems
    aztec_parameters.set("AZ_solver", "AZ_cg")
    aztec_parameters.set("AZ_output", "AZ_none")
    # settings of the ML preconditioner
    MLList = ParameterList()
    default_settings = EpetraDefaultSetter()
    default_settings.SetDefaults(MLList, "SA")
    MLList.set("max levels", 3)
    MLList.set("prec type", "MGW")
    MLList.set("smoother: type", "Chebyshev")
    MLList.set("smoother: sweeps", 2);
    # create solver
    linear_solver = MultiLevelSolver(aztec_parameters, MLList, tolerance, max_iterations)
    # only form the ML preconditioner at the first solve.
    # this speeds up solution of linear problems.
    linear_solver.SetReformPrecAtEachStep(False)

    # don't left scale the system matrix. this would destroy
    # the symmetry needed by the conjugate gradient method.
    linear_solver.SetScalingType(MLSolverScalingType.NoScaling)

    return linear_solver
