from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
from KratosMultiphysics import *
from KratosMultiphysics.TrilinosApplication import *

#
def ConstructSolver(configuration):

    import KratosMultiphysics

    if(type(configuration) != KratosMultiphysics.Parameters):
        raise Exception("input is expected to be provided as a Kratos Parameters object")

    solver_type = configuration["solver_type"].GetString()

    linear_solver = None

    #
    if(solver_type == "CGSolver" or  solver_type == "BICGSTABSolver" or solver_type == "GMRESSolver" or solver_type == "AztecSolver" ):
        linear_solver = AztecSolver(configuration)
    elif(solver_type == "MLSolver" or solver_type=="MultiLevelSolver" ):
        linear_solver = MultiLevelSolver(configuration)
    elif(solver_type == "AmgclMPISolver"):
        linear_solver = AmgclMPISolver(configuration);
    elif(solver_type == "AmgclMPISchurComplementSolver"):
        linear_solver = AmgclMPISchurComplementSolver(configuration);
    #
    elif(solver_type == "Deflated Conjugate gradient"):
        raise Exception("not implemented within trilinos")
    #
    elif(solver_type == "GMRES-UP Block"):
        raise Exception("not implemented within trilinos")
    #
    elif(solver_type == "Skyline LU factorization"):
        raise Exception("not implemented within trilinos")
    #
    elif(solver_type == "AmesosSolver"):
        linear_solver = AmesosSolver(configuration);
    #
    elif(solver_type == "SuperLUIterativeSolver"):
        raise Exception("not implemented within trilinos")
    #
    elif(solver_type == "PastixDirect"):
        raise Exception("not implemented within trilinos")
    #
    elif(solver_type == "PastixIterative"):
        raise Exception("not implemented within trilinos")
    else:
        print("*****************************************************************")
        print("Inexisting solver type. Specified::::::::::: ", solver_type)
        print(" Possibilities are:")
        print("............")
        print("*****************************************************************")
        raise Exception("please specify a correct solver")
    # else:
        # except LogicError:

    return linear_solver
