# Import Kratos
import KratosMultiphysics as KM

# Import Trilinos
import KratosMultiphysics.TrilinosApplication as KratosTrilinos

"""
@package trilinos_solver_factory

@brief Module for creating Trilinos-based linear solvers in Kratos.

This module provides essential functionality for configuring and instantiating 
    linear solvers specifically designed for the Kratos Trilinos application. 
    It ensures robustness and compatibility by:
* Checking and automatically updating deprecated solver types.
* Verifying the availability of desired solvers based on the Kratos compilation.
* Providing a utility to automatically select the most efficient direct solver.

@section Module_Functions Public and Internal Functions
* @c _CheckIfTypeIsDeprecated(config): Internal function to check and translate deprecated solver names.
* @c _CheckIfSolverIsCompiled(solver_type): Internal function to verify if the solver is available in the compilation.
* @c ConstructSolver(configuration): Primary function to construct a solver based on user configuration.
* @c CreateFastestAvailableDirectLinearSolver(): Convenience function to get the best available direct solver.

@par Example Usage
@code
import KratosMultiphysics as KM
from . import trilinos_solver_factory # Assuming this module is here

config = KM.Parameters({
"solver_type": "CGSolver",
# ... other settings
})
solver = trilinos_solver_factory.ConstructSolver(config)
@endcode
"""

def _CheckIfTypeIsDeprecated(config):
    """
    @brief Checks for deprecated solver types in the configuration.

    This function checks if the given solver type in the configuration is deprecated.
    If the type is deprecated, a warning is printed to the console, and the solver type 
        in the configuration is automatically replaced with the new, non-deprecated name 
        to ensure forward compatibility.

    @param config KM.Parameters: The configuration parameters object containing the solver settings.

    @return None

    @par Example
    @code
    _CheckIfTypeIsDeprecated(config)
    @endcode
    """
    solver_type = config["solver_type"].GetString()

    old_new_name_map = {
        "CGSolver"                      : "cg",
        "BICGSTABSolver"                : "bicgstab",
        "GMRESSolver"                   : "gmres",
        "AztecSolver"                   : "aztec",
        "MLSolver"                      : "multi_level",
        "MultiLevelSolver"              : "multi_level",
        "AmgclMPISolver"                : "amgcl",
        "AmgclMPISchurComplementSolver" : "amgcl_schur_complement",
        "AmesosSolver"                  : "amesos"
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
    @brief Checks if a specific solver type is available in the KratosTrilinos compilation.

    This function verifies whether the specified solver type, identified by a string, 
        has been compiled and is accessible within the KratosTrilinos environment.

    @param solver_type The type of the solver to check, provided as a string.

    @return @c bool: Returns @c True if the solver is both available and compiled; otherwise, @c False.

    @exception Exception Thrown if the required solver factory, which is responsible 
        for creating solvers of this type, has been completely disabled during the Kratos 
        compilation process.
    """
    # --- 1. Define Solver Groups ---
    aztec_solvers = {"aztec", "cg", "bicgstab", "gmres"}
    ml_solvers = {"multi_level"}
    
    # --- 2. Check Aztec Solvers ---
    if solver_type in aztec_solvers:
        if not hasattr(KratosTrilinos, 'AztecSolver'):
            KM.Logger.PrintWarning("Trilinos-Linear-Solver-Factory", "Trying to use Aztec-solver, which was disabled at compile time with TRILINOS_EXCLUDE_AZTEC_SOLVER!")
            return False
        return True

    # --- 3. Check MultiLevel Solvers ---
    if solver_type in ml_solvers:
        if not hasattr(KratosTrilinos, 'MultiLevelSolver'):
            KM.Logger.PrintWarning("Trilinos-Linear-Solver-Factory", "Trying to use MultiLevelSolver, which was disabled at compile time with TRILINOS_EXCLUDE_ML_SOLVER!")
            return False
        return True

    # --- 4. Check Amesos (Direct) Solvers ---
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
            KM.Logger.PrintWarning("Trilinos-Linear-Solver-Factory", "Trying to use Amesos-solver, which was disabled at compile time with TRILINOS_EXCLUDE_AMESOS_SOLVER!")
            return False

        # Case B: Only Amesos1 is compiled (Warning for performance)
        if has_amesos1 and not has_amesos2:
            KM.Logger.PrintWarning("Trilinos-Linear-Solver-Factory", "Amesos2 not compiled but Amesos is. Recommended to compile Amesos2 for better performance.")
            # 'basker' is strictly an Amesos2 solver, so it fails here
            if solver_type == "basker":
                KM.Logger.PrintWarning("Trilinos-Linear-Solver-Factory", "basker is strictly an Amesos2 solver")
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

def ConstructSolver(configuration = KM.Parameters("""{}""")):
    """
    @brief Constructs a Trilinos linear solver instance.

    This function takes a configuration object and uses the KratosTrilinos factory
        to create and return a fully configured Trilinos linear solver instance.

    @param configuration KM.Parameters: The configuration parameters object containing 
        all necessary settings for the linear solver (e.g., type, tolerance, maximum iterations).

    @return @c KratosTrilinos.TrilinosLinearSolverFactory().Create: A newly constructed 
        Trilinos linear solver object ready for use in simulations.

    @exception Exception Thrown if the input \p configuration is not a valid 
        Kratos Parameters object.

    @par Example
    @code
    solver = ConstructSolver(config)
    @endcode
    """
    if not isinstance(configuration, KM.Parameters):
        raise Exception("input is expected to be provided as a Kratos Parameters object")

    # Define solver priority list
    if not configuration.Has("solver_type"):
        linear_solvers_by_speed = KratosTrilinos.TrilinosSparseSpace.FastestDirectSolverList()
        for solver_name in linear_solvers_by_speed:
            if _CheckIfSolverIsCompiled(solver_name):
                configuration.AddEmptyValue("solver_type")
                configuration["solver_type"].SetString(solver_name)
                break

    _CheckIfTypeIsDeprecated(configuration) # for backwards-compatibility
    if not _CheckIfSolverIsCompiled(configuration["solver_type"].GetString()):
        raise Exception(configuration["solver_type"].GetString() + " solver is not available")

    return KratosTrilinos.TrilinosLinearSolverFactory().Create(configuration)

def CreateFastestAvailableDirectLinearSolver():
    """
    @brief Creates the fastest available direct linear solver.

    This function attempts to construct and return the most efficient direct linear 
        solver available in the current compilation, following a predefined priority 
        order (e.g., SuperLU, KLU, etc.) to ensure optimal performance.

    @return @c ConstructSolver: The fastest direct linear solver instance found. 
        The specific type depends on the compilation flags and availability.

    @exception Exception Thrown if, after checking all possibilities, no direct 
        linear solver could be successfully constructed due to missing libraries or 
        factory failures.

    @par Example
    @code
    fast_solver = CreateFastestAvailableDirectLinearSolver()
    @endcode
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
