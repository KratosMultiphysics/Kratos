from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
from KratosMultiphysics import *
from KratosMultiphysics.TrilinosApplication import *

#
#
#


def ConstructPreconditioner(configuration):
    preconditioner_type = "None"
    preconditioner_parameters = ParameterList()

    if hasattr(configuration, 'preconditioner'):
        preconditioner_type = configuration.preconditioner
        if(preconditioner_type == "DiagonalPreconditioner"):
            preconditioner_type = "None"
            preconditioner_parameters = ParameterList()
        elif(preconditioner_type == "ILU0"):
            preconditioner_type = "ILU"
            preconditioner_parameters = ParameterList()
            preconditioner_parameters.set("fact: drop tolerance", 1e-9)
            preconditioner_parameters.set("fact: level-of-fill", 1)
        elif(preconditioner_type == "ILUT"):
            preconditioner_type = "ILU"
            preconditioner_parameters = ParameterList()
            preconditioner_parameters.set("fact: drop tolerance", 1e-9)
            preconditioner_parameters.set("fact: level-of-fill", 1)
        #elif(preconditioner_type == "ICC"):
            #preconditioner_type = "ICC"
            #preconditioner_parameters = ParameterList()
            #preconditioner_parameters.set("fact: drop tolerance", 1e-9)
            #preconditioner_parameters.set("fact: level-of-fill", 1)
        elif(preconditioner_type == "AmesosPreconditioner"):
            preconditioner_type = "Amesos"
            preconditioner_parameters.set("amesos: solver type", "Amesos_Klu")
        else:
            raise Exception("wrong type of preconditioner")

        return [preconditioner_type, preconditioner_parameters]
    else:
        raise Exception("wrong type of preconditioner")


#
#
#
def ConstructSolver(configuration):

    ##############################################################
    ###THIS IS A VERY DIRTY HACK TO ALLOW PARAMETERS TO BE PASSED TO THE LINEAR SOLVER FACTORY
    ###TODO: clean this up!!
    if(type(configuration) == Parameters):
        solver_type = configuration["solver_type"].GetString()

        import json
        tmp = json.loads(configuration.PrettyPrintJsonString())

        class aux(object):
            pass

        configuration = aux()
        configuration.__dict__.update(tmp)
    ##############################################################

    solver_type = configuration.solver_type

    overlap_level = 0
    if hasattr(configuration, 'overlap_level'):
        overlap_level = configuration.overlap_level

    scaling = False
    if hasattr(configuration, 'scaling'):
        scaling = configuration.scaling

    verbosity = 0
    if hasattr(configuration, 'verbosity'):
        verbosity = configuration.verbosity

    linear_solver = None

    #
    if(solver_type == "Conjugate gradient"):
        [preconditioner_type, preconditioner_parameters] = ConstructPreconditioner(configuration)
        max_it = configuration.max_iteration
        tol = configuration.tolerance

        aztec_parameters = ParameterList()
        aztec_parameters.set("AZ_solver", "AZ_cg");
        aztec_parameters.set("AZ_kspace", 200);

        if(verbosity == 0):
            aztec_parameters.set("AZ_output", "AZ_none");
        else:
            aztec_parameters.set("AZ_output", verbosity);

        linear_solver = AztecSolver(aztec_parameters, preconditioner_type, preconditioner_parameters, tol, max_it, overlap_level);
    #
    elif(solver_type == "BiConjugate gradient stabilized"):
        [preconditioner_type, preconditioner_parameters] = ConstructPreconditioner(configuration)
        max_it = configuration.max_iteration
        tol = configuration.tolerance

        aztec_parameters = ParameterList()
        aztec_parameters.set("AZ_solver", "AZ_bicgstab");
        aztec_parameters.set("AZ_kspace", max_it);

        if(verbosity == 0):
            aztec_parameters.set("AZ_output", "AZ_none");
        else:
            aztec_parameters.set("AZ_output", verbosity);

        print(aztec_parameters)
        print(preconditioner_parameters)

        linear_solver = AztecSolver(aztec_parameters, preconditioner_type, preconditioner_parameters, tol, max_it, overlap_level);
    #
    elif(solver_type == "GMRES"):
        [preconditioner_type, preconditioner_parameters] = ConstructPreconditioner(configuration)
        max_it = configuration.max_iteration
        tol = configuration.tolerance

        gmres_krylov_space_dimension = 200
        if hasattr(configuration, 'gmres_krylov_space_dimension'):
            gmres_krylov_space_dimension = configuration.gmres_krylov_space_dimension

        aztec_parameters = ParameterList()
        aztec_parameters.set("AZ_solver", "AZ_gmres");
        aztec_parameters.set("AZ_kspace", gmres_krylov_space_dimension);

        if(verbosity == 0):
            aztec_parameters.set("AZ_output", "AZ_none");
        else:
            aztec_parameters.set("AZ_output", verbosity);

        linear_solver = AztecSolver(aztec_parameters, preconditioner_type, preconditioner_parameters, tol, max_it, overlap_level);


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
    elif(solver_type == "Super LU"):
        solver_parameters = ParameterList()
        linear_solver = AmesosSolver("Superludist", solver_parameters);
    #
    elif(solver_type == "SuperLUIterativeSolver"):
        raise Exception("not implemented within trilinos")
    #
    elif(solver_type == "PastixDirect"):
        raise Exception("not implemented within trilinos")
    #
    elif(solver_type == "PastixIterative"):
        raise Exception("not implemented within trilinos")
    #
    elif(solver_type == "AMGCL" or solver_type == "ML"):
        symmetric = False
        if hasattr(configuration, 'symmetric'):
            symmetric = configuration.symmetric

        max_it = configuration.max_iteration
        tol = configuration.tolerance

        max_levels = 3
        if hasattr(configuration, 'amg_max_levels'):
            amg_max_levels = configuration.amg_max_levels

        verbosity = 0
        if hasattr(configuration, 'verbosity'):
            verbosity = configuration.verbosity

        if(symmetric):
            MLList = ParameterList()
            default_settings = EpetraDefaultSetter()
            default_settings.SetDefaults(MLList, "SA");
            MLList.set("ML output", verbosity);
            MLList.set("max levels", max_levels);
            MLList.set("increasing or decreasing", "increasing");
            MLList.set("aggregation: type", "MIS");
            #MLList.set("coarse: type", "Amesos-Superludist");
            MLList.set("smoother: type", "Chebyshev");
            MLList.set("smoother: sweeps", 3);
            MLList.set("smoother: pre or post", "both");

            aztec_parameters = ParameterList()
            aztec_parameters.set("AZ_solver", "AZ_bicgstab");
            if(verbosity == 0):
                aztec_parameters.set("AZ_output", "AZ_none");
            else:
                aztec_parameters.set("AZ_output", 10);
        else:
            # settings of the ML solver
            MLList = ParameterList()

            default_settings = EpetraDefaultSetter()
            default_settings.SetDefaults(MLList, "NSSA");

            #MLList.set("ML output", 1);
            #MLList.set("coarse: max size", 10000)
            MLList.set("ML output", verbosity);
            MLList.set("max levels", max_levels);
            MLList.set("aggregation: type", "Uncoupled")
            #MLList.set("coarse: type", "Amesos-Superludist")

            aztec_parameters = ParameterList()
            aztec_parameters.set("AZ_solver", "AZ_gmres");
            if(verbosity == 0):
                aztec_parameters.set("AZ_output", "AZ_none");
            else:
                aztec_parameters.set("AZ_output", verbosity);

        linear_solver = MultiLevelSolver(aztec_parameters, MLList, tol, max_it);
        return linear_solver
    else:
        print("*****************************************************************")
        print("Inexisting solver type. Possibilities are:")
        print("............")
        print("*****************************************************************")
        err
    # else:
        # except LogicError:

    if(scaling == False):
        return linear_solver
    else:
        linear_solver.SetScalingType(AztecScalingType.LeftScaling)
        return linear_solver
