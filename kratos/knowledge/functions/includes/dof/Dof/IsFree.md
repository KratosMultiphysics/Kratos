---
type: C++ Method
title: IsFree
resource: includes/dof.h#L378-L381
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/future/solving_strategies/builders/elimination_builder/EliminationBuilder/AllocateLinearSystemConstraints
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/future/solving_strategies/schemes/implicit_scheme/ImplicitScheme/UpdateConstraintsOnlyDofs
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/future/solving_strategies/schemes/static_scheme/StaticScheme/Update
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/solving_strategies/builder_and_solvers/p_multigrid/augmented_lagrange_constraint_assembler/ApplyDirichletConditions
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/solving_strategies/builder_and_solvers/residualbased_elimination_builder_and_solver_with_constraints/ResidualBasedEliminationBuilderAndSolverWithConstraints/AssembleRHSWithoutConstraints
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/solving_strategies/convergencecriterias/displacement_criteria/DisplacementCriteria/IsFreeAndLocalDof
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/solving_strategies/convergencecriterias/mixed_generic_criteria/MixedGenericCriteria/GetNormValues
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/solving_strategies/convergencecriterias/residual_criteria/ResidualCriteria/IsFreeAndLocalDof
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/solving_strategies/schemes/residual_based_adjoint_static_scheme/ResidualBasedAdjointStaticScheme/Initialize
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/solving_strategies/schemes/residualbased_incremental_aitken_static_scheme/ResidualBasedIncrementalAitkenStaticScheme/Update
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/sources/test_dof/KRATOS_TEST_CASE_IN_SUITE
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/sources/test_node/KRATOS_TEST_CASE_IN_SUITE
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/dof_updater/DofUpdater/UpdateDofs
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/dof_updater/DofUpdater/AssignDofs
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`bool IsFree() const noexcept`

# Called by

- [AllocateLinearSystemConstraints](../../../../functions/future/solving_strategies/builders/elimination_builder/EliminationBuilder/AllocateLinearSystemConstraints.md)
- [UpdateConstraintsOnlyDofs](../../../../functions/future/solving_strategies/schemes/implicit_scheme/ImplicitScheme/UpdateConstraintsOnlyDofs.md)
- [Update](../../../../functions/future/solving_strategies/schemes/static_scheme/StaticScheme/Update.md)
- [ApplyDirichletConditions](../../../../functions/solving_strategies/builder_and_solvers/p_multigrid/augmented_lagrange_constraint_assembler/ApplyDirichletConditions.md)
- [AssembleRHSWithoutConstraints](../../../../functions/solving_strategies/builder_and_solvers/residualbased_elimination_builder_and_solver_with_constraints/ResidualBasedEliminationBuilderAndSolverWithConstraints/AssembleRHSWithoutConstraints.md)
- [IsFreeAndLocalDof](../../../../functions/solving_strategies/convergencecriterias/displacement_criteria/DisplacementCriteria/IsFreeAndLocalDof.md)
- [GetNormValues](../../../../functions/solving_strategies/convergencecriterias/mixed_generic_criteria/MixedGenericCriteria/GetNormValues.md)
- [IsFreeAndLocalDof](../../../../functions/solving_strategies/convergencecriterias/residual_criteria/ResidualCriteria/IsFreeAndLocalDof.md)
- [Initialize](../../../../functions/solving_strategies/schemes/residual_based_adjoint_static_scheme/ResidualBasedAdjointStaticScheme/Initialize.md)
- [Update](../../../../functions/solving_strategies/schemes/residualbased_incremental_aitken_static_scheme/ResidualBasedIncrementalAitkenStaticScheme/Update.md)
- [KRATOS_TEST_CASE_IN_SUITE](../../../../functions/tests/cpp_tests/sources/test_dof/KRATOS_TEST_CASE_IN_SUITE.md)
- [KRATOS_TEST_CASE_IN_SUITE](../../../../functions/tests/cpp_tests/sources/test_node/KRATOS_TEST_CASE_IN_SUITE.md)
- [UpdateDofs](../../../../functions/utilities/dof_updater/DofUpdater/UpdateDofs.md)
- [AssignDofs](../../../../functions/utilities/dof_updater/DofUpdater/AssignDofs.md)