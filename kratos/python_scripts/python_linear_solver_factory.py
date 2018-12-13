from __future__ import print_function, absolute_import, division # makes KratosMultiphysics backward compatible with python 2.6 and 2.7

def ConstructSolver(configuration):
    import KratosMultiphysics as KM

    if(type(configuration) != KM.Parameters):
        raise Exception("input is expected to be provided as a Kratos Parameters object")

    solver_type = configuration["solver_type"].GetString()

    # dict specifying in which applications the linear-solvers are defined
    # it is necessary to import the applications in which the linear-solvers
    # are defined, because otherwise they are not registered in the kernel
    # and hence cannot be created by the C-linear-solver-factory
    linear_solver_apps = {
        "ExternalSolversApplication" : [
            "GMRESSolver",
            "SuperLUSolver",
            "Super LU",
            "Super_LU",
            "SuperLUIterativeSolver",
            "PastixSolver",
            "complex_pastix_solver"
        ],
        "EigenSolversApplication" : [
            "eigen_sparse_lu",
            "eigen_pardiso_llt",
            "eigen_pardiso_ldlt",
            "eigen_pardiso_lu",
            "complex_eigen_sparse_lu",
            "complex_eigen_pardiso_lu"
        ]
    }

    for app_name, linear_solve_names in linear_solver_apps.items():
        if solver_type in linear_solve_names:
            # import the Application in which the linear solver is defined
            __import__("KratosMultiphysics." + app_name)
            break

    if KM.ComplexLinearSolverFactory().Has(solver_type):
        return KM.ComplexLinearSolverFactory().Create(configuration)
    else:
        return KM.LinearSolverFactory().Create(configuration)
