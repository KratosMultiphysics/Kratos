from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7

import KratosMultiphysics as KM
import KratosMultiphysics.TrilinosApplication as KratosTrilinos

def _CheckIfTypeIsDeprecated(config):
    '''function to translate old/deprecated names to new names
    needed for backwards-compatibility
    '''
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
    if solver_type in ["aztec", "cg", "bicgstab", "gmres"] and not hasattr(KratosTrilinos, 'AztecSolver'):
        raise Exception('Trying to use Aztec-solver, which was disabled at compile time with "TRILINOS_EXCLUDE_AZTEC_SOLVER"!')
    if solver_type in ["amesos", "klu", "super_lu_dist", "mumps"] and not hasattr(KratosTrilinos, 'AmesosSolver'):
        raise Exception('Trying to use Amesos-solver, which was disabled at compile time with "TRILINOS_EXCLUDE_AMESOS_SOLVER"!')
    if solver_type in ["multi_level"] and not hasattr(KratosTrilinos, 'MultiLevelSolver'):
        raise Exception('Trying to use MultiLevelSolver-solver, which was diasbled at compile time with "TRILINOS_EXCLUDE_ML_SOLVER"!')

def ConstructSolver(configuration):
    if not isinstance(configuration, KM.Parameters):
        raise Exception("input is expected to be provided as a Kratos Parameters object")

    _CheckIfTypeIsDeprecated(configuration) # for backwards-compatibility
    _CheckIfSolverIsCompiled(configuration["solver_type"].GetString())

    return KratosTrilinos.TrilinosLinearSolverFactory().Create(configuration)
