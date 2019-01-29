from __future__ import print_function, absolute_import, division # makes KratosMultiphysics backward compatible with python 2.6 and 2.7

import KratosMultiphysics as KM

def __DeprecatedApplicationImport(solver_type):
    # dict specifying in which applications the linear-solvers are defined
    # it is necessary to import the applications in which the linear-solvers
    # are defined, because otherwise they are not registered in the kernel
    # and hence cannot be created by the C-linear-solver-factory
    # NOTE: this is only for backwards-compatibility!!!
    # the correct way is to specify the application in which the linear solver
    # is defined, e.g. "solver_type" : "ExternalSolversApplication.SuperLUSolver"

    linear_solver_apps = {
        "ExternalSolversApplication" : [
            "GMRESSolver",
            "SuperLUSolver",
            "Super_LU",
            "SuperLUIterativeSolver",
            "PastixSolver",
            "PastixComplexSolver",
            "FEASTSolver"
        ],
        "EigenSolversApplication" : [
            "SparseLUSolver",
            "eigen_sparse_lu",
            "SparseQRSolver",
            "ComplexSparseLUSolver",
            "ComplexSparseQRSolver",
            "PardisoLLTSolver",
            "eigen_pardiso_llt",
            "ComplexPardisoLLTSolver",
            "PardisoLDLTSolver",
            "eigen_pardiso_ldlt",
            "eigen_pardiso_lu",
            "ComplexPardisoLDLTSolver",
            "complex_eigen_sparse_lu",
            "complex_eigen_pardiso_lu"
            "ComplexPardisoLUSolver"
            "PardisoLUSolver"
        ]
    }

    for app_name, linear_solver_names in linear_solver_apps.items():
        if solver_type in linear_solver_names:
            depr_msg  = 'The linear-solver "' + solver_type + '" is defined in the "' + app_name +  '"\n'
            depr_msg += 'Please specify the "solver_type" including the name of the application:\n'
            depr_msg += '"' + app_name + '.' + solver_type + '"'
            depr_msg += '\nPlease update your settings accordingly, the current settings are deprecated!'
            KM.Logger.PrintWarning('DEPRECATION-WARNING', depr_msg)

            from KratosMultiphysics import kratos_utilities as kratos_utils
            if not kratos_utils.IsApplicationAvailable(app_name):
                err_msg  = 'Trying to use the linear-solver "' + solver_type
                err_msg += '"\nThis solver is defined in the "' + app_name
                err_msg += '" which is not compiled'
                raise Exception(err_msg)
            # import the Application in which the linear solver is defined
            __import__("KratosMultiphysics." + app_name)
            break


def ConstructSolver(configuration):
    if(type(configuration) != KM.Parameters):
        raise Exception("input is expected to be provided as a Kratos Parameters object")

    solver_type = configuration["solver_type"].GetString()

    # remove unused "KratosMultiphysics.
    if solver_type.startswith("KratosMultiphysics."):
        solver_type = solver_type[19:]

    if "Application." in solver_type: # the module in which the solver is implemented was specified
        splitted_name = solver_type.split(".")
        if len(splitted_name) != 2:
            raise NameError('The "solver_type" has to consist in "ApplicationName.SolverType"')
        app_name = splitted_name[0]
        solver_type = splitted_name[1]
        configuration["solver_type"].SetString(solver_type)
        __import__("KratosMultiphysics." + app_name)
    else:
        __DeprecatedApplicationImport(solver_type)

    if KM.ComplexLinearSolverFactory().Has(solver_type):
        KM.Logger.PrintInfo("Linear-Solver-Factory",\
            "Constructing a complex linear-solver")
        return KM.ComplexLinearSolverFactory().Create(configuration)
    else:
        KM.Logger.PrintInfo("Linear-Solver-Factory",\
            "Constructing a regular (non-complex) linear-solver")
        return KM.LinearSolverFactory().Create(configuration)
