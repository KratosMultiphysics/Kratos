---
type: C++ Function
title: CalculateEntityAuxiliaryVariableContributions
resource: solving_strategies/schemes/residual_based_adjoint_bossak_scheme.h#L903-L927
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/utilities/openmp_utils/ThisThread
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/solving_strategies/schemes/residual_based_adjoint_bossak_scheme/ResidualBasedAdjointBossakScheme/CheckAndResizeThreadStorage
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/includes/amatrix_interface/prod
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/solving_strategies/schemes/residual_based_adjoint_bossak_scheme/ResidualBasedAdjointBossakScheme/CalculateAuxiliaryVariableContributions
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void CalculateEntityAuxiliaryVariableContributions( TEntityType& rCurrentEntity, LocalSystemVectorType& rAdjointAuxiliaryValues, const ProcessInfo& rProcessInfo)`

# Calls

- [ThisThread](../../../../functions/utilities/openmp_utils/ThisThread.md)
- [CheckAndResizeThreadStorage](../../../../functions/solving_strategies/schemes/residual_based_adjoint_bossak_scheme/ResidualBasedAdjointBossakScheme/CheckAndResizeThreadStorage.md)
- [prod](../../../../functions/includes/amatrix_interface/prod.md)

# Called by

- [CalculateAuxiliaryVariableContributions](../../../../functions/solving_strategies/schemes/residual_based_adjoint_bossak_scheme/ResidualBasedAdjointBossakScheme/CalculateAuxiliaryVariableContributions.md)