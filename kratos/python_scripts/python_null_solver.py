class PythonNullSolver:
    """
    Null-object implementation of the solver interface.

    This class is used to allow running an AnalysisStage without an actual solver.
    Instead of checking everywhere whether a solver exists, `_GetSolver()` always
    returns a valid object. When no real solver is provided, this `_NullSolver`
    instance is used.

    Behavior
    --------
    - Any public method call (e.g. `Initialize`, `SolveSolutionStep`, etc.)
      is accepted and silently ignored (no-op).
    - This ensures that the execution pipeline remains unchanged even when
      no solver is present.

    Important notes
    ---------------
    - Only *public* methods are allowed. Access to private/protected attributes
      (names starting with "_") will raise an AttributeError.
    - This class intentionally hides errors such as calling non-existing methods
      (e.g. due to typos). This is acceptable because the solver API is stable.
    - This class must remain lightweight and stateless.

    Implementation details
    ----------------------
    - `__getattr__` dynamically returns a no-op function for any undefined method.
    - `__slots__ = ()` prevents dynamic attribute creation and keeps memory usage minimal.
    """

    __slots__ = ()

    def __getattr__(self, name):
        if name.startswith("_"):
            raise AttributeError(f"PythonNullSolver: access to '{name}' is not allowed")
        return lambda *a, **k: None

# Singleton instance used across all AnalysisStage objects
# This enforces statelessness and avoids unnecessary repeated allocations
PYTHON_NULL_SOLVER = PythonNullSolver()