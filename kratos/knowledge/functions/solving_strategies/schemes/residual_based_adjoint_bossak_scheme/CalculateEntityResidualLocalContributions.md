---
type: C++ Function
title: CalculateEntityResidualLocalContributions
resource: solving_strategies/schemes/residual_based_adjoint_bossak_scheme.h#L678-L689
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/utilities/openmp_utils/ThisThread
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/includes/amatrix_interface/prod
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/solving_strategies/schemes/residual_based_adjoint_bossak_scheme/ResidualBasedAdjointBossakScheme/CalculateResidualLocalContributions
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void CalculateEntityResidualLocalContributions( TEntityType& rCurrentEntity, LocalSystemMatrixType& rLHS_Contribution, LocalSystemVectorType& rRHS_Contribution, const ProcessInfo& rCurrentProcessInfo)`

# Calls

- [ThisThread](../../../../functions/utilities/openmp_utils/ThisThread.md)
- [prod](../../../../functions/includes/amatrix_interface/prod.md)

# Called by

- [CalculateResidualLocalContributions](../../../../functions/solving_strategies/schemes/residual_based_adjoint_bossak_scheme/ResidualBasedAdjointBossakScheme/CalculateResidualLocalContributions.md)