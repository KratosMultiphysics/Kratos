---
type: C++ Method
title: Assemble_AdjointVars
resource: solving_strategies/schemes/residual_based_adjoint_bossak_scheme.h#L1309-L1336
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/solving_strategies/schemes/residual_based_adjoint_bossak_scheme/ResidualBasedAdjointBossakScheme/UpdateTimeSchemeAdjoints
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/solving_strategies/schemes/residual_based_adjoint_bossak_scheme/ResidualBasedAdjointBossakScheme/UpdateAuxiliaryVariable
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`static void Assemble_AdjointVars(const std::vector<const VariableData*>& rVariables, Communicator& rComm)`

# Called by

- [UpdateTimeSchemeAdjoints](../../../../../functions/solving_strategies/schemes/residual_based_adjoint_bossak_scheme/ResidualBasedAdjointBossakScheme/UpdateTimeSchemeAdjoints.md)
- [UpdateAuxiliaryVariable](../../../../../functions/solving_strategies/schemes/residual_based_adjoint_bossak_scheme/ResidualBasedAdjointBossakScheme/UpdateAuxiliaryVariable.md)