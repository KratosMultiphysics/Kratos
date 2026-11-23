---
type: C++ Method
title: UpdateTimeSchemeAdjoints
resource: solving_strategies/schemes/residual_based_adjoint_bossak_scheme.h#L946-L996
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/solving_strategies/schemes/residual_based_adjoint_bossak_scheme/ResidualBasedAdjointBossakScheme/GatherVariables
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/solving_strategies/schemes/residual_based_adjoint_bossak_scheme/ResidualBasedAdjointBossakScheme/SetToZero_AdjointVars
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/solving_strategies/schemes/residual_based_adjoint_bossak_scheme/UpdateEntityTimeSchemeContributions
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/solving_strategies/schemes/residual_based_adjoint_bossak_scheme/ResidualBasedAdjointBossakScheme/Assemble_AdjointVars
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/solving_strategies/schemes/residual_based_adjoint_bossak_scheme/ResidualBasedAdjointBossakScheme/Update
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void UpdateTimeSchemeAdjoints(ModelPart& rModelPart)`

# Calls

- [GatherVariables](../../../../../functions/solving_strategies/schemes/residual_based_adjoint_bossak_scheme/ResidualBasedAdjointBossakScheme/GatherVariables.md)
- [SetToZero_AdjointVars](../../../../../functions/solving_strategies/schemes/residual_based_adjoint_bossak_scheme/ResidualBasedAdjointBossakScheme/SetToZero_AdjointVars.md)
- [UpdateEntityTimeSchemeContributions](../../../../../functions/solving_strategies/schemes/residual_based_adjoint_bossak_scheme/UpdateEntityTimeSchemeContributions.md)
- [Assemble_AdjointVars](../../../../../functions/solving_strategies/schemes/residual_based_adjoint_bossak_scheme/ResidualBasedAdjointBossakScheme/Assemble_AdjointVars.md)

# Called by

- [Update](../../../../../functions/solving_strategies/schemes/residual_based_adjoint_bossak_scheme/ResidualBasedAdjointBossakScheme/Update.md)