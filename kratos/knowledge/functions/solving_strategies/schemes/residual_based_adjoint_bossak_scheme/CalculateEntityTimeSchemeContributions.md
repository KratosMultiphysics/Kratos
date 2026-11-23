---
type: C++ Function
title: CalculateEntityTimeSchemeContributions
resource: solving_strategies/schemes/residual_based_adjoint_bossak_scheme.h#L853-L888
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
  - target: functions/solving_strategies/schemes/residual_based_adjoint_bossak_scheme/ResidualBasedAdjointBossakScheme/CalculateTimeSchemeContributions
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void CalculateEntityTimeSchemeContributions( TEntityType& rCurrentEntity, LocalSystemVectorType& rAdjointTimeSchemeValues2, LocalSystemVectorType& rAdjointTimeSchemeValues3, const ProcessInfo& rProcessInfo)`

# Calls

- [ThisThread](../../../../functions/utilities/openmp_utils/ThisThread.md)
- [CheckAndResizeThreadStorage](../../../../functions/solving_strategies/schemes/residual_based_adjoint_bossak_scheme/ResidualBasedAdjointBossakScheme/CheckAndResizeThreadStorage.md)
- [prod](../../../../functions/includes/amatrix_interface/prod.md)

# Called by

- [CalculateTimeSchemeContributions](../../../../functions/solving_strategies/schemes/residual_based_adjoint_bossak_scheme/ResidualBasedAdjointBossakScheme/CalculateTimeSchemeContributions.md)