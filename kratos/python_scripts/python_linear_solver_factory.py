from __future__ import print_function, absolute_import, division # makes KratosMultiphysics backward compatible with python 2.6 and 2.7

import KratosMultiphysics as KM
from KratosMultiphysics import kratos_utilities as kratos_utils
from importlib import import_module

def __DeprecatedApplicationImport(solver_type):
    # dict specifying in which applications the linear-solvers are defined
    # it is necessary to import the applications in which the linear-solvers
    # are defined, because otherwise they are not registered in the kernel
    # and hence cannot be created by the C-linear-solver-factory
    # NOTE: this is only for backwards-compatibility!!!
    # the correct way is to specify the application in which the linear solver
    # is defined, e.g. "solver_type" : "ExternalSolversApplication.super_lu"

    linear_solver_apps = {
        "ExternalSolversApplication" : [
            "gmres",
            "super_lu",
            "super_lu_iterative",
            "pastix",
            "pastix_complex",
            "feast"
        ],
        "EigenSolversApplication" : [
            "sparse_lu",
            "sparse_lu_complex",
            "sparse_qr",
            "sparse_qr_complex",
            "pardiso_llt",
            "pardiso_ldlt",
            "pardiso_lu",
            "pardiso_lu_complex"
        ]
    }

    for app_name, linear_solver_names in linear_solver_apps.items():
        if solver_type in linear_solver_names:
            depr_msg  = 'DEPRECATION-WARNING:\nThe linear-solver "' + solver_type + '" is defined in the "' + app_name +  '"\n'
            depr_msg += 'Please specify the "solver_type" including the name of the application:\n'
            depr_msg += '"' + app_name + '.' + solver_type + '"'
            depr_msg += '\nPlease update your settings accordingly, the current settings are deprecated!'
            KM.Logger.PrintWarning('Linear-Solver-Factory', depr_msg)

            if not kratos_utils.CheckIfApplicationsAvailable(app_name):
                err_msg  = 'Trying to use the linear-solver "' + solver_type
                err_msg += '"\nThis solver is defined in the "' + app_name
                err_msg += '" which is not compiled'
                raise Exception(err_msg)
            # import the Application in which the linear solver is defined
            import_module("KratosMultiphysics." + app_name)
            break

def __CheckIfSolverTypeIsDeprecated(config):
    '''function to translate old/deprecated solver-names to new names
    needed for backwards-compatibility
    '''

    solver_type = config["solver_type"].GetString()
    splitted_solver_type = solver_type.split(".") # in case the name of the solver is "EigenSolversApp.SparseLUSolver"

    real_solver_type = splitted_solver_type[-1]

    # solvers from Core
    old_new_name_map = {
        "CGSolver"                     : "cg",
        "BICGSTABSolver"               : "bicgstab",
        "DeflatedCGSolver"             : "deflated_cg",
        "TFQMRSolver"                  : "tfqmr",
        "SkylineLUFactorizationSolver" : "skyline_lu_factorization",
        "AMGCL"                        : "amgcl",
        "AMGCLSolver"                  : "amgcl",
        "AMGCL_NS_Solver"              : "amgcl_ns",
        "ScalingSolver"                : "scaling",
        "SkylineLUComplexSolver"       : "skyline_lu_complex",
        "complex_skyline_lu_solver"    : "skyline_lu_complex"
    }

    # solvers from ExternalSolversApp
    old_new_name_map.update({
        "GMRESSolver"            : "gmres",
        "Super_LU"               : "super_lu",
        "SuperLUSolver"          : "super_lu",
        "SuperLUIterativeSolver" : "super_lu_iterative",
        "PastixSolver"           : "pastix",
        "PastixComplexSolver"    : "pastix_complex",
        "FEASTSolver"            : "feast"
    })

    # solvers from EigenSolversApp
    old_new_name_map.update({
        "SparseLUSolver"           : "sparse_lu",
        "eigen_sparse_lu"          : "sparse_lu",
        "ComplexSparseLUSolver"    : "sparse_lu_complex",
        "complex_eigen_sparse_lu"  : "sparse_lu_complex",
        "SparseQRSolver"           : "sparse_qr",
        "ComplexSparseQRSolver"    : "sparse_qr_complex",
        "PardisoLLTSolver"         : "pardiso_llt",
        "eigen_pardiso_llt"        : "pardiso_llt",
        "PardisoLDLTSolver"        : "pardiso_ldlt",
        "eigen_pardiso_ldlt"       : "pardiso_ldlt",
        "PardisoLUSolver"          : "pardiso_lu",
        "eigen_pardiso_lu"         : "pardiso_lu",
        "ComplexPardisoLUSolver"   : "pardiso_lu_complex",
        "complex_eigen_pardiso_lu" : "pardiso_lu_complex"
    })

    if real_solver_type in old_new_name_map:
        new_name = old_new_name_map[real_solver_type]
        depr_msg  = 'DEPRECATION-WARNING: \nUsing a deprecated "solver_type"!\n'
        depr_msg += 'Replace "' + real_solver_type + '" with "' + new_name + '"'
        KM.Logger.PrintWarning("Linear-Solver-Factory", depr_msg)
        splitted_solver_type[-1] = new_name
        config["solver_type"].SetString(".".join(splitted_solver_type))

def __CheckIfPreconditionerTypeIsDeprecated(config):
    '''function to translate old/deprecated preconditioner-names to new names
    needed for backwards-compatibility
    '''

    if config.Has("preconditioner_type"):
        preconditioner_type = config["preconditioner_type"].GetString()

        old_new_name_map = {
            "None"                   : "none",
            "DiagonalPreconditioner" : "diagonal",
            "ILU0Preconditioner"     : "ilu0",
            "ILUPreconditioner"      : "ilu"
        }

        if preconditioner_type in old_new_name_map:
            new_name = old_new_name_map[preconditioner_type]
            depr_msg  = 'DEPRECATION-WARNING: \nUsing a deprecated "preconditioner_type"!\n'
            depr_msg += 'Replace "' + preconditioner_type + '" with "' + new_name + '"'
            KM.Logger.PrintWarning("Linear-Solver-Factory", depr_msg)
            config["preconditioner_type"].SetString(new_name)

def ConstructSolver(configuration):
    if(type(configuration) != KM.Parameters):
        raise Exception("input is expected to be provided as a Kratos Parameters object")

    __CheckIfSolverTypeIsDeprecated(configuration)
    __CheckIfPreconditionerTypeIsDeprecated(configuration)

    solver_type = configuration["solver_type"].GetString()

    # remove unused "KratosMultiphysics.
    if solver_type.startswith("KratosMultiphysics."):
        solver_type = solver_type[19:]

    if "Application." in solver_type: # the module in which the solver is implemented was specified
        splitted_name = solver_type.split(".")
        if len(splitted_name) != 2:
            raise NameError('The "solver_type" has to consist in "ApplicationName.solver_type"')
        app_name = splitted_name[0]
        # the following is only needed for the check in the ComplexLinearSolverFactory
        # note that the solver-configuration is NOT modified
        solver_type = splitted_name[1]
        import_module("KratosMultiphysics." + app_name)
    else:
        __DeprecatedApplicationImport(solver_type)

    if KM.ComplexLinearSolverFactory().Has(solver_type):
        KM.Logger.PrintInfo("Linear-Solver-Factory", "Constructing a complex linear-solver")
        return KM.ComplexLinearSolverFactory().Create(configuration)
    else:
        KM.Logger.PrintInfo("Linear-Solver-Factory", "Constructing a regular (non-complex) linear-solver")
        return KM.LinearSolverFactory().Create(configuration)

def CreateFastestAvailableDirectLinearSolver():
    # Using a default linear solver (selecting the fastest one available)
    if kratos_utils.CheckIfApplicationsAvailable("EigenSolversApplication"):
        from KratosMultiphysics import EigenSolversApplication
    elif kratos_utils.CheckIfApplicationsAvailable("ExternalSolversApplication"):
        from KratosMultiphysics import ExternalSolversApplication

    linear_solvers_by_speed = [
        "pardiso_lu", # EigenSolversApplication (if compiled with Intel-support)
        "sparse_lu",  # EigenSolversApplication
        "pastix",     # ExternalSolversApplication (if Pastix is included in compilation)
        "super_lu",   # ExternalSolversApplication
        "skyline_lu_factorization" # in Core, always available, but slow
    ]

    for solver_name in linear_solvers_by_speed:
        if KM.LinearSolverFactory().Has(solver_name):
            linear_solver_configuration = KM.Parameters("""{ "solver_type" : "%s"}""" % solver_name)

            KM.Logger.PrintInfo('Linear-Solver-Factory', 'Creating "{}" as fastest available direct solver'.format(solver_name))
            return KM.LinearSolverFactory().Create(linear_solver_configuration)

    raise Exception("Linear-Solver could not be constructed!")
