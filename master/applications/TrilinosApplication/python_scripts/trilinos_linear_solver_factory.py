# Import Kratos
import KratosMultiphysics as KM

# Import Trilinos
import KratosMultiphysics.TrilinosApplication as KratosTrilinos

"""
This module provides functionality for creating a linear solver for use with the Trilinos application.
It includes checking and converting deprecated solver types, checking if the desired solver is compiled, and creating the fastest available direct linear solver.

Importing necessary Kratos modules:
- `KratosMultiphysics` as `KM`
- `KratosMultiphysics.TrilinosApplication` as `KratosTrilinos`

Functions:
- `_CheckIfTypeIsDeprecated(config)`: Checks and translates deprecated solver names to new names for backward compatibility.
- `_CheckIfSolverIsCompiled(solver_type)`: Checks if the desired solver is compiled.
- `ConstructSolver(configuration)`: Constructs the solver based on the given configuration.
- `CreateFastestAvailableDirectLinearSolver()`: Creates the fastest available direct linear solver.

Example:
    config = KM.Parameters({
        "solver_type": "CGSolver"
    })
    solver = ConstructSolver(config)
"""

def _CheckIfTypeIsDeprecated(config):
    """
    This function checks if the given solver type in the configuration is deprecated.
    If the type is deprecated, a warning is printed and the solver type in the configuration is replaced with the new name.

    Args:
        config (KM.Parameters): The configuration parameters for the solver.

    Example:
        _CheckIfTypeIsDeprecated(config)
    """
    solver_type = config["solver_type"].GetString()

    old_new_name_map = {
        "CGSolver" : "cg",
        "BICGSTABSolver" : "bicgstab",
        "GMRESSolver" : "gmres",
        "AztecSolver" : "aztec",
        "MLSolver" : "multi_level",
        "MultiLevelSolver" : "multi_level",
        "AmgclMPISolver" : "amgcl",
        "AmgclMPISchurComplementSolver" : "amgcl_schur_complement",
        "AmesosSolver" : "amesos"
    }

    if solver_type in old_new_name_map:
        new_name = old_new_name_map[solver_type]
        depr_msg  = '\nDEPRECATION-WARNING: using a deprecated "solver_type"!\n'
        depr_msg += 'Replace "' + solver_type + '" with "' + new_name + '"'
        import KratosMultiphysics as KM
        KM.Logger.PrintWarning("Trilinos-Linear-Solver-Factory", depr_msg)
        config["solver_type"].SetString(new_name)

def _CheckIfSolverIsCompiled(solver_type):
    """
    This function checks if the specified solver type is compiled.

    Args:
        solver_type (str): The type of the solver to check.

    Raises:
        Exception: If the solver type is not compiled.

    Example:
        _CheckIfSolverIsCompiled("aztec")
    """
    if solver_type in ["aztec", "cg", "bicgstab", "gmres"] and not hasattr(KratosTrilinos, 'AztecSolver'):
        raise Exception('Trying to use Aztec-solver, which was disabled at compile time with "TRILINOS_EXCLUDE_AZTEC_SOLVER"!')
    if solver_type in ["amesos", "klu", "super_lu_dist", "mumps"]:
        has_amesos = hasattr(KratosTrilinos, 'AmesosSolver')
        has_amesos_2 = hasattr(KratosTrilinos, 'Amesos2Solver')
        if not has_amesos_2:
            if has_amesos:
                KM.Logger.PrintWarning("Trilinos-Linear-Solver-Factory", "Amesos2 not compiled but Amesos it is, recommended for better performance")
            else:
                raise Exception('Trying to use Amesos-solver, which was disabled at compile time with "TRILINOS_EXCLUDE_AMESOS_SOLVER"!')
    if solver_type in ["multi_level"] and not hasattr(KratosTrilinos, 'MultiLevelSolver'):
        raise Exception('Trying to use MultiLevelSolver-solver, which was diasbled at compile time with "TRILINOS_EXCLUDE_ML_SOLVER"!')

def ConstructSolver(configuration):
    """
    Constructs and returns a Trilinos linear solver based on the given configuration.

    Args:
        configuration (KM.Parameters): The configuration parameters for the solver.

    Returns:
        KratosTrilinos.TrilinosLinearSolverFactory().Create: A Trilinos linear solver.

    Raises:
        Exception: If the input is not a Kratos Parameters object.

    Example:
        solver = ConstructSolver(config)
    """
    if not isinstance(configuration, KM.Parameters):
        raise Exception("input is expected to be provided as a Kratos Parameters object")

    _CheckIfTypeIsDeprecated(configuration) # for backwards-compatibility
    _CheckIfSolverIsCompiled(configuration["solver_type"].GetString())

    return KratosTrilinos.TrilinosLinearSolverFactory().Create(configuration)

def CreateFastestAvailableDirectLinearSolver():
    """
    Creates and returns the fastest available direct linear solver, based on a predefined order and availability.

    TODO: A proper study of speed must be done, particularly depending of system size would be interesting

    In Trilinos, the speed of direct solvers like MUMPS, SuperLU_DIST, KLU, and Basker can depend on various factors including the size and sparsity of the matrix, the architecture of the computer system, and the specific characteristics of the problem being solved. Below are some general observations about these solvers:

    1. MUMPS (MUltifrontal Massively Parallel Sparse Direct Solver):
    - Parallel solver, handles large problems well.
    - Can be used on distributed-memory machines.
    - May have higher memory requirements.

    2. SuperLU_DIST:
    - Parallel solver, also designed for distributed-memory machines.
    - Often used for large, sparse linear systems.
    - Has good performance on a wide range of problems.

    3. KLU:
    - A serial solver, typically used for smaller problems.
    - Works well for circuit simulation problems and other problems with a similar structure.
    - Generally not suitable for large distributed-memory parallel computations.

    4. Basker:
    - Also a serial solver, suitable for smaller problems.
    - May not perform as well on larger problems or on distributed-memory systems.

    Here are some considerations for choosing a solver:

    - If you are working on a distributed-memory parallel machine and dealing with large problems, you might want to consider MUMPS or SuperLU_DIST.
    - For smaller problems or problems with structure similar to circuit simulation problems, KLU might be a good choice.
    - Basker might be suitable for smaller problems where other solvers are not performing well.

    Returns:
        ConstructSolver: The fastest available direct linear solver.

    Raises:
        Exception: If no linear solver could be constructed.

    Example:
        fast_solver = CreateFastestAvailableDirectLinearSolver()
    """
    linear_solvers_by_speed = [
        "mumps", # TODO: Which first?, MUMPS of SuperLUDist?
        "super_lu_dist", 
        "klu",
        "basker"
    ]

    registered_names_solvers_amesos = {
        "mumps"         : "Amesos_Mumps",
        "super_lu_dist" : "Amesos_Superludist", 
        "klu"           : "Amesos_Klu"
    }

    registered_names_solvers_amesos2 = {
        "mumps"         : "amesos2_mumps",
        "super_lu_dist" : "amesos2_superludist", 
        "klu"           : "amesos2_klu2",
        "basker"        : "basker"
    }

    # Settings
    settings = KM.Parameters("""{"solver_type" : ""}""")
    for solver_name in linear_solvers_by_speed:
        if hasattr(KratosTrilinos, 'Amesos2Solver'):
            registered_name = registered_names_solvers_amesos2[solver_name]
            if KratosTrilinos.Amesos2Solver.HasSolver(registered_name):
                settings["solver_type"].SetString("amesos2")
                settings.AddEmptyValue("amesos2_solver_type")
                settings["amesos2_solver_type"].SetString(registered_name)
                return ConstructSolver(settings)
        elif hasattr(KratosTrilinos, 'AmesosSolver'):
            registered_name = registered_names_solvers_amesos[solver_name]
            if KratosTrilinos.AmesosSolver.HasSolver(registered_name):
                settings["solver_type"].SetString(solver_name)
                return ConstructSolver(settings)

    raise Exception("Linear-Solver could not be constructed!")
