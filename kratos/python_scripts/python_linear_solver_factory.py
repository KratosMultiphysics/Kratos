import KratosMultiphysics as KM
from KratosMultiphysics import kratos_utilities as kratos_utils
from importlib import import_module

def ConstructSolver(Settings):
    """Construct a linear solver using the factory specified in ``Settings``.

    Reads the ``kratos_module`` and ``factory_module`` entries from
    ``Settings`` to locate and instantiate the appropriate
    ``PythonLinearSolverFactory`` subclass, then delegates to its
    ``CreateLinearSolver`` method.

    Missing factory-related keys (``kratos_module``, ``factory_module``,
    ``solver_type``) are filled in with default values pointing to
    ``PythonLinearSolverFactory``.

    Args:
        Settings (KratosMultiphysics.Parameters): Parameters object that must
            contain at least a ``"solver_type"`` entry. The optional
            ``"kratos_module"`` and ``"factory_module"`` entries control which
            factory class is used.

    Returns:
        LinearSolver: A Kratos linear solver instance.

    Raises:
        Exception: If ``Settings`` is not a
            ``KratosMultiphysics.Parameters`` object.
    """
    # Check input
    if type(Settings) != KM.Parameters:
        raise Exception("input is expected to be provided as a Kratos Parameters object")

    # Default factory settings (do NOT inject these into Settings to avoid polluting the solver config)
    kratos_module_name = "KratosMultiphysics"
    if Settings.Has("kratos_module"):
        kratos_module_name = Settings["kratos_module"].GetString()
        Settings.RemoveValue("kratos_module")
    factory_module_name = "python_linear_solver_factory.PythonLinearSolverFactory"
    if Settings.Has("factory_module"):
        factory_module_name = Settings["factory_module"].GetString()
        Settings.RemoveValue("factory_module")

    # Generate the python module for the factory
    if not kratos_module_name.startswith("KratosMultiphysics"):
        kratos_module_name = "KratosMultiphysics." + kratos_module_name # ensure the correct prefix is added if not provided
    module_path, class_name = factory_module_name.rsplit(".", 1)
    factory_module = import_module(kratos_module_name + "." + module_path)

    # Generate the factory
    factory_prototype = getattr(factory_module, class_name)
    factory_instance = factory_prototype()

    # Generate the solver from the factory and return it
    return factory_instance.CreateLinearSolver(Settings)

def CreateFastestAvailableDirectLinearSolver(Settings = KM.Parameters("{}")):
    """Create the fastest available direct linear solver using the configured factory.

    Reads the ``kratos_module`` and ``factory_module`` entries from
    ``Settings`` (defaulting to ``PythonLinearSolverFactory``) to locate and
    instantiate the appropriate factory class, then delegates to its
    ``CreateFastestAvailableDirectLinearSolver`` method.

    Args:
        Settings (KratosMultiphysics.Parameters, optional): Parameters object
            that may contain ``"kratos_module"`` and ``"factory_module"``
            entries to override the default factory. Defaults to an empty
            Parameters object, which selects ``PythonLinearSolverFactory``.

    Returns:
        LinearSolver: The fastest available direct Kratos linear solver.

    Raises:
        Exception: If no suitable direct solver can be found.
    """
    # Fill in missing factory keys from defaults (do NOT inject into Settings)
    kratos_module_name = Settings["kratos_module"].GetString() if Settings.Has("kratos_module") else "KratosMultiphysics"
    factory_module_name = Settings["factory_module"].GetString() if Settings.Has("factory_module") else "python_linear_solver_factory.PythonLinearSolverFactory"

    # Resolve and import the factory module
    if not kratos_module_name.startswith("KratosMultiphysics"):
        kratos_module_name = "KratosMultiphysics." + kratos_module_name
    module_path, class_name = factory_module_name.rsplit(".", 1)
    factory_module = import_module(kratos_module_name + "." + module_path)

    # Instantiate the factory and delegate
    factory_prototype = getattr(factory_module, class_name)
    factory_instance = factory_prototype()
    return factory_instance.CreateFastestAvailableDirectLinearSolver()

class PythonLinearSolverFactory(object):
    """Factory class for creating linear solvers in Kratos.

    Provides methods to construct linear solvers (both regular and complex)
    from a Parameters configuration object, and to select the fastest available
    direct solver on the current system.
    """
    def CreateLinearSolver(self, configuration):
        """Create a linear solver from a Parameters configuration object.

        Automatically determines whether to construct a complex or regular
        (non-complex) linear solver based on the registered solver types.
        If the solver belongs to an external application, that application
        is imported first.

        Args:
            configuration (KratosMultiphysics.Parameters): Parameters object
                containing at least a ``"solver_type"`` entry. Optionally may
                contain an ``"inner_solver_settings"`` block for preconditioners
                or nested solvers.

        Returns:
            LinearSolver: A Kratos linear solver instance.

        Raises:
            Exception: If ``configuration`` is not a
                ``KratosMultiphysics.Parameters`` object.
        """
        if type(configuration) != KM.Parameters:
            raise Exception("input is expected to be provided as a Kratos Parameters object")

        solver_type = self.__GetSolverTypeAndImportApplication(configuration["solver_type"].GetString())

        if configuration.Has("inner_solver_settings"):
            inner_solver_type = configuration["inner_solver_settings"]["solver_type"].GetString()
            self.__GetSolverTypeAndImportApplication(inner_solver_type)

        if KM.ComplexLinearSolverFactory().Has(solver_type):
            KM.Logger.PrintInfo(self.__GetLabel(), "Constructing a complex linear-solver")
            return KM.ComplexLinearSolverFactory().Create(configuration)
        else:
            KM.Logger.PrintInfo(self.__GetLabel(), "Constructing a regular (non-complex) linear-solver")
            return KM.LinearSolverFactory().Create(configuration)

    def CreateFastestAvailableDirectLinearSolver(self):
        """Create the fastest available direct linear solver on the current system.

        Iterates over the list of direct solvers ordered by speed
        (``UblasSparseSpace.FastestDirectSolverList``) and returns the first
        one that is registered in the ``LinearSolverFactory``. If the
        ``LinearSolversApplication`` is available it is loaded beforehand so
        that its solvers are registered.

        Returns:
            LinearSolver: The fastest available direct Kratos linear solver.

        Raises:
            Exception: If no suitable direct solver can be found.
        """
        # Using a default linear solver (selecting the fastest one available)
        if kratos_utils.CheckIfApplicationsAvailable("LinearSolversApplication"):
            from KratosMultiphysics import LinearSolversApplication

        linear_solvers_by_speed = KM.UblasSparseSpace.FastestDirectSolverList()

        for solver_name in linear_solvers_by_speed:
            if KM.LinearSolverFactory().Has(solver_name):
                linear_solver_configuration = KM.Parameters("""{ "solver_type" : "%s"}""" % solver_name)

                KM.Logger.PrintInfo(self.__GetLabel(), 'Creating "{}" as fastest available direct solver'.format(solver_name))
                return KM.LinearSolverFactory().Create(linear_solver_configuration)

        raise Exception("Linear-Solver could not be constructed!")

    def __GetLabel(self):
        """Return the label used for logger messages.

        Returns:
            str: The label string ``"Linear-Solver-Factory"``.
        """
        return "Linear-Solver-Factory"

    def __GetSolverTypeAndImportApplication(self, solver_type):
        """Resolve the solver type string and import its application if needed.

        Strips the leading ``"KratosMultiphysics."`` prefix when present. If
        the solver type includes an application qualifier
        (``"ApplicationName.solver_type"``), the corresponding Kratos
        application module is imported so that its solvers are registered in
        the factory.

        Args:
            solver_type (str): Raw solver type string, e.g.
                ``"LinearSolversApplication.sparse_lu"`` or
                ``"KratosMultiphysics.LinearSolversApplication.sparse_lu"``.

        Returns:
            str: The bare solver type name without application or namespace
            prefix, e.g. ``"sparse_lu"``.

        Raises:
            NameError: If the qualified name does not follow the expected
                ``"ApplicationName.solver_type"`` format.
        """
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