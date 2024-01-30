import KratosMultiphysics as KM
from KratosMultiphysics import kratos_utilities as kratos_utils
from importlib import import_module


def ConstructSolver(configuration):
    if(type(configuration) != KM.Parameters):
        raise Exception("input is expected to be provided as a Kratos Parameters object")

    solver_type = __GetSolverTypeAndImportApplication(configuration["solver_type"].GetString())

    if configuration.Has("inner_solver_settings"):
        inner_solver_type = configuration["inner_solver_settings"]["solver_type"].GetString()
        __GetSolverTypeAndImportApplication(inner_solver_type)

    if KM.ComplexLinearSolverFactory().Has(solver_type):
        KM.Logger.PrintInfo("Linear-Solver-Factory", "Constructing a complex linear-solver")
        return KM.ComplexLinearSolverFactory().Create(configuration)
    else:
        KM.Logger.PrintInfo("Linear-Solver-Factory", "Constructing a regular (non-complex) linear-solver")
        return KM.LinearSolverFactory().Create(configuration)


def CreateFastestAvailableDirectLinearSolver():
    # Using a default linear solver (selecting the fastest one available)
    if kratos_utils.CheckIfApplicationsAvailable("LinearSolversApplication"):
        from KratosMultiphysics import LinearSolversApplication

    linear_solvers_by_speed = [
        "pardiso_lu", # LinearSolversApplication (if compiled with Intel-support)
        "sparse_lu",  # LinearSolversApplication
        "skyline_lu_factorization" # in Core, always available, but slow
    ]

    for solver_name in linear_solvers_by_speed:
        if KM.LinearSolverFactory().Has(solver_name):
            linear_solver_configuration = KM.Parameters("""{ "solver_type" : "%s"}""" % solver_name)

            KM.Logger.PrintInfo('Linear-Solver-Factory', 'Creating "{}" as fastest available direct solver'.format(solver_name))
            return KM.LinearSolverFactory().Create(linear_solver_configuration)

    raise Exception("Linear-Solver could not be constructed!")


def __GetSolverTypeAndImportApplication(solver_type):
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

    return solver_type

