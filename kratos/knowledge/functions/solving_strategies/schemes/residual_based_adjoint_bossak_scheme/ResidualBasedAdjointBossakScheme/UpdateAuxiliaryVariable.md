---
type: C++ Method
title: UpdateAuxiliaryVariable
resource: solving_strategies/schemes/residual_based_adjoint_bossak_scheme.h#L1093-L1113
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
  - target: functions/solving_strategies/schemes/residual_based_adjoint_bossak_scheme/UpdateEntityAuxiliaryVariableContributions
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/solving_strategies/schemes/residual_based_adjoint_bossak_scheme/ResidualBasedAdjointBossakScheme/Assemble_AdjointVars
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/solving_strategies/schemes/residual_based_adjoint_bossak_scheme/ResidualBasedAdjointBossakScheme/FinalizeSolutionStep
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void UpdateAuxiliaryVariable(ModelPart& rModelPart)`

# Calls

- [GatherVariables](../../../../../functions/solving_strategies/schemes/residual_based_adjoint_bossak_scheme/ResidualBasedAdjointBossakScheme/GatherVariables.md)
- [SetToZero_AdjointVars](../../../../../functions/solving_strategies/schemes/residual_based_adjoint_bossak_scheme/ResidualBasedAdjointBossakScheme/SetToZero_AdjointVars.md)
- [UpdateEntityAuxiliaryVariableContributions](../../../../../functions/solving_strategies/schemes/residual_based_adjoint_bossak_scheme/UpdateEntityAuxiliaryVariableContributions.md)
- [Assemble_AdjointVars](../../../../../functions/solving_strategies/schemes/residual_based_adjoint_bossak_scheme/ResidualBasedAdjointBossakScheme/Assemble_AdjointVars.md)

# Called by

- [FinalizeSolutionStep](../../../../../functions/solving_strategies/schemes/residual_based_adjoint_bossak_scheme/ResidualBasedAdjointBossakScheme/FinalizeSolutionStep.md)