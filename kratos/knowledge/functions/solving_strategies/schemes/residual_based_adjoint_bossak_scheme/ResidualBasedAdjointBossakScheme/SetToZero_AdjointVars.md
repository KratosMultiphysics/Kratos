---
type: C++ Method
title: SetToZero_AdjointVars
resource: solving_strategies/schemes/residual_based_adjoint_bossak_scheme.h#L1280-L1307
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/utilities/variable_utils/SetHistoricalVariableToZero
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/solving_strategies/schemes/residual_based_adjoint_bossak_scheme/ResidualBasedAdjointBossakScheme/UpdateTimeSchemeAdjoints
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/solving_strategies/schemes/residual_based_adjoint_bossak_scheme/ResidualBasedAdjointBossakScheme/UpdateAuxiliaryVariable
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`static void SetToZero_AdjointVars(const std::vector<const VariableData*>& rVariables, ModelPart::NodesContainerType& rNodes)`

# Calls

- [SetHistoricalVariableToZero](../../../../../functions/utilities/variable_utils/SetHistoricalVariableToZero.md)

# Called by

- [UpdateTimeSchemeAdjoints](../../../../../functions/solving_strategies/schemes/residual_based_adjoint_bossak_scheme/ResidualBasedAdjointBossakScheme/UpdateTimeSchemeAdjoints.md)
- [UpdateAuxiliaryVariable](../../../../../functions/solving_strategies/schemes/residual_based_adjoint_bossak_scheme/ResidualBasedAdjointBossakScheme/UpdateAuxiliaryVariable.md)