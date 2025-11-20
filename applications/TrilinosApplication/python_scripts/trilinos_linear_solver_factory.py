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
    Checks if the specified solver type is compiled in KratosTrilinos.
    
    Args:
        solver_type (str): The type of the solver to check.
        
    Returns:
        bool: True if the solver is available and compiled.
        
    Raises:
        Exception: If the required solver factory is completely disabled at compile time.
    """
    
    # --- 1. Define Solver Groups ---
    aztec_solvers = {"aztec", "cg", "bicgstab", "gmres"}
    ml_solvers = {"multi_level"}
    
    # --- 2. Check Aztec Solvers ---
    if solver_type in aztec_solvers:
        if not hasattr(KratosTrilinos, 'AztecSolver'):
            raise Exception('Trying to use Aztec-solver, which was disabled at compile time with "TRILINOS_EXCLUDE_AZTEC_SOLVER"!')
        return True

    # --- 3. Check MultiLevel Solvers ---
    if solver_type in ml_solvers:
        if not hasattr(KratosTrilinos, 'MultiLevelSolver'):
            raise Exception('Trying to use MultiLevelSolver, which was disabled at compile time with "TRILINOS_EXCLUDE_ML_SOLVER"!')
        return True

    # --- 4. Check Amesos (Direct) Solvers ---
    # This set includes generic names and specific implementation names
    amesos_solvers = {
        "amesos", "amesos2", "klu", "klu2", 
        "super_lu_dist", "super_lu_dist2", 
        "mumps", "mumps2", "basker"
    }

    if solver_type in amesos_solvers:
        has_amesos1 = hasattr(KratosTrilinos, 'AmesosSolver')
        has_amesos2 = hasattr(KratosTrilinos, 'Amesos2Solver')

        # Case A: Neither is compiled
        if not has_amesos1 and not has_amesos2:
            raise Exception('Trying to use Amesos-solver, which was disabled at compile time with "TRILINOS_EXCLUDE_AMESOS_SOLVER"!')

        # Case B: Only Amesos1 is compiled (Warning for performance)
        if has_amesos1 and not has_amesos2:
            KM.Logger.PrintWarning("Trilinos-Linear-Solver-Factory", "Amesos2 not compiled but Amesos is. Recommended to compile Amesos2 for better performance.")
            # 'basker' is strictly an Amesos2 solver, so it fails here
            if solver_type == "basker":
                return False

        # Case C: specific internal solver check
        # Map inputs to (required_factory, internal_trilinos_name)
        requires_amesos2 = solver_type.endswith("2") or solver_type == "basker"
        
        # Check availability based on required factory
        if requires_amesos2:
            if not has_amesos2: 
                return False
            
            # Map external name -> internal Trillinos name
            map_amesos2 = {
                "mumps2": "amesos2_mumps",
                "super_lu_dist2": "amesos2_superludist",
                "klu2": "amesos2_klu2",
                "basker": "basker"
            }
            # If it's a generic name like "amesos2", we assume True, otherwise check specific availability
            internal_name = map_amesos2.get(solver_type)
            if internal_name and not KratosTrilinos.Amesos2Solver.HasSolver(internal_name):
                return False
        else: # requires Amesos1
            if not has_amesos1:
                return False
                
            map_amesos1 = {
                "mumps": "Amesos_Mumps",
                "super_lu_dist": "Amesos_Superludist",
                "klu": "Amesos_Klu"
            }
            internal_name = map_amesos1.get(solver_type)
            if internal_name and not KratosTrilinos.AmesosSolver.HasSolver(internal_name):
                return False

    return True

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

    Returns:
        ConstructSolver: The fastest available direct linear solver.

    Raises:
        Exception: If no linear solver could be constructed.

    Example:
        fast_solver = CreateFastestAvailableDirectLinearSolver()
    """
    linear_solvers_by_speed = KratosTrilinos.TrilinosSparseSpace.FastestDirectSolverList()

    # Settings
    settings = KM.Parameters("""{"solver_type" : ""}""")
    for solver_name in linear_solvers_by_speed:
        if _CheckIfSolverIsCompiled(solver_name):
            KM.Logger.PrintInfo("Trilinos-Linear-Solver-Factory", 'Using "' + solver_name + '" as the fastest available direct linear solver.')
            settings["solver_type"].SetString(solver_name)
            return ConstructSolver(settings)

    raise Exception("Linear-Solver could not be constructed!")
