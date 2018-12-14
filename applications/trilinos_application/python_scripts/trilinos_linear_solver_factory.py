from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7

def _CheckIfTypeIsDeprecated(config):
    '''function to translate old/deprecated names to new names
    needed for backwards-compatibility
    '''
    solver_type = config["solver_type"].GetString()

    old_new_name_map = {
        "CGSolver" : "aztec_cg",
        "BICGSTABSolver" : "aztec_bicgstab",
        "GMRESSolver" : "aztec_gmres",
        "AztecSolver" : "aztec",
        "MLSolver" : "multi_level",
        "MultiLevelSolver" : "multi_level",
        "AmgclMPISolver" : "amgcl_mpi",
        "AmgclMPISchurComplementSolver" : "amgcl_mpi_schur_complement",
        "AmesosSolver" : "amesos"
    }

    if solver_type in old_new_name_map:
        new_name = old_new_name_map[solver_type]
        depr_msg  = 'WARNING, using a deprecated "solver_type"!\n'
        depr_msg += 'Replace "' + solver_type + '" with "' + new_name + '"'
        import KratosMultiphysics as KM
        KM.Logger.PrintWarning("Trilinos-Linear-Solver", depr_msg)
        config["solver_type"].SetString(new_name)

def ConstructSolver(configuration):

    import KratosMultiphysics as KM
    import KratosMultiphysics.TrilinosApplication as KratosTrilinos

    if(type(configuration) != KM.Parameters):
        raise Exception("input is expected to be provided as a Kratos Parameters object")

    _CheckIfTypeIsDeprecated(configuration) # for backwards-compatibility

    return KratosTrilinos.TrilinosLinearSolverFactory().Create(configuration)
