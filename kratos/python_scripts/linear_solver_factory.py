from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7

from KratosMultiphysics import *

def ConstructPreconditioner(configuration):
    if hasattr(configuration, 'preconditioner_type'):
        preconditioner_type = configuration.preconditioner_type
        if(preconditioner_type == "None"):
            return None
        else:
            if(preconditioner_type == "DiagonalPreconditioner"):
                return DiagonalPreconditioner()
            elif(preconditioner_type == "ILU0Preconditioner"):
                return ILU0Preconditioner()
            else:
                print("Preconditioner type not found. Returning None")
                return None
    else:
        return None


#
#
#
def ConstructSolver(configuration):

    depr_msg  = '"kratos/python_scripts/linear_solver_factory.py" is deprecated and will be removed\n'
    depr_msg += 'Please use "kratos/python_scripts/python_linear_solver_factory.py" instead!'
    Logger.PrintWarning('DEPRECATION-WARNING', depr_msg)

    params = 0
    ##############################################################
    ###THIS IS A VERY DIRTY HACK TO ALLOW PARAMETERS TO BE PASSED TO THE LINEAR SOLVER FACTORY
    ###TODO: clean this up!!
    if(type(configuration) == Parameters):
        from KratosMultiphysics import python_linear_solver_factory
        return python_linear_solver_factory.ConstructSolver(configuration)
        #solver_type = configuration["solver_type"].GetString()

        #import json
        #tmp = json.loads(configuration.PrettyPrintJsonString())

        #params = configuration

        #class aux(object):
            #pass

        #configuration = aux()
        #configuration.__dict__.update(tmp)
    ##############################################################






    solver_type = configuration.solver_type

    scaling = False
    if hasattr(configuration, 'scaling'):
        scaling = configuration.scaling

    linear_solver = None

    #
    if(solver_type == "Conjugate gradient" or solver_type == "Conjugate_gradient" or solver_type == "CGSolver"):
        precond = ConstructPreconditioner(configuration)
        max_it = configuration.max_iteration
        tol = configuration.tolerance
        if(precond is None):
            linear_solver = CGSolver(tol, max_it)
        else:
            linear_solver = CGSolver(tol, max_it, precond)
    #
    elif(solver_type == "BiConjugate gradient stabilized" or solver_type == "BiConjugate_gradient_stabilized" or solver_type == "BICGSTABSolver"):
        precond = ConstructPreconditioner(configuration)
        max_it = configuration.max_iteration
        tol = configuration.tolerance
        if(precond is None):
            linear_solver = BICGSTABSolver(tol, max_it)
        else:
            linear_solver = BICGSTABSolver(tol, max_it, precond)
    #
    elif(solver_type == "GMRES" or solver_type == "GMRESSolver"):
        import KratosMultiphysics.ExternalSolversApplication
        precond = ConstructPreconditioner(configuration)
        max_it = configuration.max_iteration
        tol = configuration.tolerance
        if(precond is None):
            linear_solver = KratosMultiphysics.ExternalSolversApplication.GMRESSolver(
                tol, max_it)
        else:
            linear_solver = KratosMultiphysics.ExternalSolversApplication.GMRESSolver(
                tol, max_it, precond)
    #
    elif(solver_type == "Deflated Conjugate gradient" or solver_type == "Deflated_Conjugate_gradient" or solver_type == "DeflatedCGSolver"):
        max_it = configuration.max_iteration
        tol = configuration.tolerance

        assume_constant_structure = False
        if hasattr(configuration, 'assume_constant_structure'):
            assume_constant_structure = configuration.assume_constant_structure

        max_reduced_size = 1000
        if hasattr(configuration, 'assume_constant_structure'):
            max_reduced_size = configuration.max_reduced_size

        linear_solver = DeflatedCGSolver(
            tol,
            max_it,
            assume_constant_structure,
            max_reduced_size)
    #
    elif(solver_type == "Skyline LU factorization" or solver_type == "Skyline_LU_factorization" or solver_type == "SkylineLUFactorizationSolver"):
        linear_solver = SkylineLUFactorizationSolver()
    #
    elif(solver_type == "Super LU" or solver_type == "Super_LU" or solver_type == "SuperLUSolver"):
        import KratosMultiphysics.ExternalSolversApplication
        linear_solver = KratosMultiphysics.ExternalSolversApplication.SuperLUSolver()
    #
    elif(solver_type == "SuperLUIterativeSolver"):
        import KratosMultiphysics.ExternalSolversApplication
        tol = configuration.tolerance
        max_it = configuration.max_iteration

        if(hasattr(configuration, "gmres_krylov_space_dimension")):
            restart = configuration.gmres_krylov_space_dimension
        else:
            print("WARNING: restart not specitifed, setting it to the number of iterations")
            restart = max_it

        if(hasattr(configuration, "DropTol")):
            DropTol = configuration.DropTol
        else:
            print("WARNING: DropTol not specified, setting it to 1e-4")
            DropTol = 1e-4

        if(hasattr(configuration, "FillTol")):
            FillTol = configuration.FillTol
        else:
            print("WARNING: FillTol not specified, setting it to 1e-2")
            FillTol = 1e-2

        if(hasattr(configuration, "ilu_level_of_fill")):
            ilu_level_of_fill = configuration.ilu_level_of_fill
        else:
            print("WARNING: level of fill not specified, setting it to 10")
            ilu_level_of_fill = 10
        linear_solver = KratosMultiphysics.ExternalSolversApplication.SuperLUIterativeSolver(
            tol, max_it, restart, DropTol, FillTol, ilu_level_of_fill)

    #
    elif(solver_type == "PastixDirect" or solver_type == "PastixSolver"):
        import KratosMultiphysics.ExternalSolversApplication
        is_symmetric = False
        if hasattr(configuration, 'is_symmetric'):
            configuration.is_symmetric
        verbosity = 0
        if hasattr(configuration, 'verbosity'):
            verbosity = configuration.verbosity
        linear_solver = KratosMultiphysics.ExternalSolversApplication.PastixSolver(
            verbosity, is_symmetric)
    #
    elif(solver_type == "PastixIterative"):
        import KratosMultiphysics.ExternalSolversApplication
        tol = configuration.tolerance
        max_it = configuration.max_iteration
        restart = configuration.gmres_krylov_space_dimension
        ilu_level_of_fill = configuration.ilu_level_of_fill
        is_symmetric = configuration.is_symmetric
        verbosity = 0
        if hasattr(configuration, 'verbosity'):
            verbosity = configuration.verbosity
        linear_solver = KratosMultiphysics.ExternalSolversApplication.PastixSolver(
            tol, restart, ilu_level_of_fill, verbosity, is_symmetric)
    #
    elif(solver_type == "AMGCL"):

        if(params == 0): #old style construction
            if hasattr(configuration, 'preconditioner_type'):
                if(configuration.preconditioner_type != "None"):
                    print("WARNING: preconditioner specified in preconditioner_type will not be used as it is not compatible with the AMGCL solver")

            max_it = configuration.max_iteration
            tol = configuration.tolerance

            verbosity = 0
            if hasattr(configuration, 'verbosity'):
                verbosity = configuration.verbosity

            if hasattr(configuration, 'smoother_type'):
                smoother_type = configuration.smoother_type  # options are DAMPED_JACOBI, ILU0, SPAI
                if(smoother_type == "ILU0"):
                    amgcl_smoother = AMGCLSmoother.ILU0
                elif(smoother_type == "DAMPED_JACOBI"):
                    amgcl_smoother = AMGCLSmoother.DAMPED_JACOBI
                elif(smoother_type == "SPAI0"):
                    amgcl_smoother = AMGCLSmoother.SPAI0
                elif(smoother_type == "GAUSS_SEIDEL"):
                    amgcl_smoother = AMGCLSmoother.GAUSS_SEIDEL
                else:
                    print("ERROR: smoother_type shall be one of \"ILU0\", \"DAMPED_JACOBI\", \"SPAI0\", got \"{0}\".\n\"ILU0\" will be used.".format(smoother_type))
                    amgcl_smoother = AMGCLSmoother.ILU0
            else:
                print("WARNING: smoother_type not prescribed for AMGCL solver, setting it to ILU0")
                amgcl_smoother = AMGCLSmoother.ILU0

            if hasattr(configuration, 'krylov_type'):
                krylov_type = configuration.krylov_type
                if(krylov_type == "GMRES"):
                    amgcl_krylov_type = AMGCLIterativeSolverType.GMRES
                if(krylov_type == "LGMRES"):
                    amgcl_krylov_type = AMGCLIterativeSolverType.LGMRES
                elif(krylov_type == "BICGSTAB"):
                    amgcl_krylov_type = AMGCLIterativeSolverType.BICGSTAB
                elif(krylov_type == "CG"):
                    amgcl_krylov_type = AMGCLIterativeSolverType.CG
                elif(krylov_type == "BICGSTAB2"):
                    amgcl_krylov_type = AMGCLIterativeSolverType.BICGSTAB2
                elif(krylov_type == "BICGSTAB_WITH_GMRES_FALLBACK"):
                    amgcl_krylov_type = AMGCLIterativeSolverType.BICGSTAB_WITH_GMRES_FALLBACK
                else:
                    print("ERROR: krylov_type shall be one of \"GMRES\", \"LGMRES\",\"BICGSTAB\", \"CG\", \"BICGSTAB2\", \"BICGSTAB_WITH_GMRES_FALLBACK\", got \"{0}\".\n\"GMRES\" will be used".format(krylov_type))
                    amgcl_krylov_type = AMGCLIterativeSolverType.GMRES
            else:
                print("WARNING: krylov_type not prescribed for AMGCL solver, setting it to GMRES")
                amgcl_krylov_type = AMGCLIterativeSolverType.GMRES



            if hasattr(configuration, 'coarsening_type'):
                coarsening_type = configuration.coarsening_type
                if(coarsening_type == "RUGE_STUBEN"):
                    amgcl_coarsening_type = AMGCLCoarseningType.RUGE_STUBEN
                elif(coarsening_type == "AGGREGATION"):
                    amgcl_coarsening_type = AMGCLCoarseningType.AGGREGATION
                elif(coarsening_type == "SA"):
                    amgcl_coarsening_type = AMGCLCoarseningType.SA
                elif(coarsening_type == "SA_EMIN"):
                    amgcl_coarsening_type = AMGCLCoarseningType.SA_EMIN
                else:
                    print("ERROR: coarsening_type shall be one of \"RUGE_STUBEN\", \"AGGREGATION\", \"SA\", \"SA_EMIN\", got \"{0}\".\n\"AGGREGATION\" will be used".format(krylov_type))
                    amgcl_coarsening_type = AMGCLCoarseningType.AGGREGATION
            else:
                amgcl_coarsening_type = AMGCLCoarseningType.AGGREGATION

            if hasattr(configuration, 'gmres_krylov_space_dimension'):
                m = configuration.gmres_krylov_space_dimension
            else:
                m = max_it


            if hasattr(configuration, 'provide_coordinates'):
                provide_coordinates = configuration.provide_coordinates
            else:
                provide_coordinates = False

            linear_solver = AMGCLSolver(
                amgcl_smoother, amgcl_krylov_type, amgcl_coarsening_type, tol, max_it, verbosity, m, provide_coordinates)

        else: ##construction by parameters
            linear_solver = AMGCLSolver(params)

    elif(solver_type == "AMGCL_NS_Solver"):
        if(params == 0): #old style construction
            params = Parameters("""
                {
                                       "solver_type" : "AMGCL_NS_Solver",
                                       "krylov_type" : "gmres",
                                       "velocity_block_preconditioner" :
                                        {
                                            "krylov_type" : "gmres",
                                            "tolerance" : 1e-3,
                                            "preconditioner_type" : "ilu0"
                                        },
                                        "pressure_block_preconditioner" :
                                        {
                                            "krylov_type" : "bicgstab",
                                            "tolerance" : 1e-2,
                                            "preconditioner_type" : "spai0"
                                        },
                                       "tolerance" : 1e-6,
                                       "gmres_krylov_space_dimension": 50,
                                       "coarsening_type": "aggregation",
                                       "max_iteration": 50,
                                       "verbosity" : 1,
                                       "scaling": ,
                                       "coarse_enough" : 5000
                                   }
                """)
        scaling = params["scaling"].GetBool()

        linear_solver = AMGCL_NS_Solver(params)

    elif (solver_type == "Parallel MKL Pardiso" or solver_type == "Parallel_MKL_Pardiso"):
        # emulating the solvers of the MKLSolversApplication through the EigenSolversApplication
        Logger.PrintWarning("LinearSolverFactor", "Solver Parallel_MKL_Pardiso is deprecated,\
        please use it through the EigenSolversApplication (see the Readme in the Application)")
        import EigenSolversApplication
        params = Parameters("""{}""")
        linear_solver = EigenSolversApplication.PardisoLUSolver(params)

    else:
        print("*****************************************************************")
        print("Inexisting solver type. Possibilities are:")
        print("Conjugate gradient")
        print("BiConjugate gradient stabilized")
        print("GMRES")
        print("Deflated Conjugate gradient")
        print("AMGCL")
        print("GMRES-UP Block")
        print("Skyline LU factorization")
        print("Super LU (requires ExternalSolversApplication)")
        print("SuperLUIterativeSolver (requires ExternalSolversApplication)")
        print("PastixDirect (requires ExternalSolversApplication + shall be habilitated at compilation time)")
        print("PastixIterative (requires ExternalSolversApplication + shall be habilitated at compilation time)")
        print("Parallel MKL Pardiso (requires EigenSolversApplication with MKL enabled)")
        print("*****************************************************************")
        raise RuntimeError(" Wrong Solver Definition ")
    # else:
        # except LogicError:

    if(scaling == False):
        return linear_solver
    else:
        return ScalingSolver(linear_solver, True)
