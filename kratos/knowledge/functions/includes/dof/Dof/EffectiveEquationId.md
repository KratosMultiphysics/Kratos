---
type: C++ Method
title: EffectiveEquationId
resource: includes/dof.h#L313-L316
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/future/solving_strategies/builders/block_builder/BlockBuilder/ApplyBlockBuildDirichletConditions
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/future/solving_strategies/builders/builder/Builder/SetUpMasterSlaveConstraintsGraph
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/future/solving_strategies/builders/builder/Builder/CalculateSolutionVector
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/future/solving_strategies/builders/elimination_builder/EliminationBuilder/AllocateLinearSystemConstraints
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/future/solving_strategies/schemes/implicit_scheme/ImplicitScheme/BuildMasterSlaveConstraints
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/future/solving_strategies/schemes/implicit_scheme/ImplicitScheme/UpdateConstraintsOnlyDofs
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/sources/test_dof/KRATOS_TEST_CASE_IN_SUITE
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/utilities/test_dof_array_utilities/KRATOS_TEST_CASE_IN_SUITE
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`EquationIdType EffectiveEquationId() const noexcept`

# Called by

- [ApplyBlockBuildDirichletConditions](../../../../functions/future/solving_strategies/builders/block_builder/BlockBuilder/ApplyBlockBuildDirichletConditions.md)
- [SetUpMasterSlaveConstraintsGraph](../../../../functions/future/solving_strategies/builders/builder/Builder/SetUpMasterSlaveConstraintsGraph.md)
- [CalculateSolutionVector](../../../../functions/future/solving_strategies/builders/builder/Builder/CalculateSolutionVector.md)
- [AllocateLinearSystemConstraints](../../../../functions/future/solving_strategies/builders/elimination_builder/EliminationBuilder/AllocateLinearSystemConstraints.md)
- [BuildMasterSlaveConstraints](../../../../functions/future/solving_strategies/schemes/implicit_scheme/ImplicitScheme/BuildMasterSlaveConstraints.md)
- [UpdateConstraintsOnlyDofs](../../../../functions/future/solving_strategies/schemes/implicit_scheme/ImplicitScheme/UpdateConstraintsOnlyDofs.md)
- [KRATOS_TEST_CASE_IN_SUITE](../../../../functions/tests/cpp_tests/sources/test_dof/KRATOS_TEST_CASE_IN_SUITE.md)
- [KRATOS_TEST_CASE_IN_SUITE](../../../../functions/tests/cpp_tests/utilities/test_dof_array_utilities/KRATOS_TEST_CASE_IN_SUITE.md)