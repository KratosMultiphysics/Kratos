from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
from KratosMultiphysics import *
from KratosMultiphysics.TrilinosApplication import *

def LinearSolver(tolerance, max_iterations):
    # settings for the iterative solver
    aztec_parameters = ParameterList()
    aztec_parameters.set("AZ_solver", "AZ_gmres")
    # aztec_parameters.set("AZ_output","AZ_none");
    aztec_parameters.set("AZ_output", 10)

    # settings of the ML solver
    MLList = ParameterList()

    MultiLevelSolver.SetDefaults(MLList, "NSSA")

    MLList.set("ML output", 1)
    MLList.set("coarse: max size", 10000)
    MLList.set("max levels", 3)
    # MLList.set("increasing or decreasing","increasing");
    # MLList.set("null space: add default vectors", 1);
    MLList.set("aggregation: type", "Uncoupled")
    # MLList.set("smoother: sweeps",3);
    # MLList.set("smoother: pre or post", "both")
    # MLList.set("ML output",0);
    MLList.set("coarse: type", "Amesos-Superludist")
    # MLList.set("smoother: ifpack type", "ILU");
    # MLList.set("smoother: ifpack overlap", 0);
    # MLList.SetSublistIntValue("smoother: ifpack list","fact: level-of-fill", 5);
    # MLList.set("coarse: sweeps", 3)
    # MLList.set("coarse: pre or post", "both")

    # MLList.set("print unused",1)
    # MLList.setboolvalue("energy minimization: enable",0)
    # MLList.set("aggregation: damping factor",0.0)

    linear_solver = MultiLevelSolver(aztec_parameters, MLList, tolerance, max_iterations)

    return linear_solver
