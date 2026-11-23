---
type: C++ Method
title: CheckAndResizeThreadStorage
resource: solving_strategies/schemes/residual_based_adjoint_bossak_scheme.h#L612-L645
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/utilities/openmp_utils/ThisThread
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/solving_strategies/schemes/residual_based_adjoint_bossak_scheme/ResidualBasedAdjointBossakScheme/CalculateSystemContributions
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/solving_strategies/schemes/residual_based_adjoint_bossak_scheme/CalculateEntityTimeSchemeContributions
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/solving_strategies/schemes/residual_based_adjoint_bossak_scheme/CalculateEntityAuxiliaryVariableContributions
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`virtual void CheckAndResizeThreadStorage(unsigned SystemSize)`

# Calls

- [ThisThread](../../../../../functions/utilities/openmp_utils/ThisThread.md)

# Called by

- [CalculateSystemContributions](../../../../../functions/solving_strategies/schemes/residual_based_adjoint_bossak_scheme/ResidualBasedAdjointBossakScheme/CalculateSystemContributions.md)
- [CalculateEntityTimeSchemeContributions](../../../../../functions/solving_strategies/schemes/residual_based_adjoint_bossak_scheme/CalculateEntityTimeSchemeContributions.md)
- [CalculateEntityAuxiliaryVariableContributions](../../../../../functions/solving_strategies/schemes/residual_based_adjoint_bossak_scheme/CalculateEntityAuxiliaryVariableContributions.md)