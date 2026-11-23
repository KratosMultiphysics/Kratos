---
type: C++ Method
title: GatherVariables
resource: solving_strategies/schemes/residual_based_adjoint_bossak_scheme.h#L1258-L1278
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/utilities/openmp_utils/ThisThread
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/solving_strategies/schemes/residual_based_adjoint_bossak_scheme/ResidualBasedAdjointBossakScheme/Check
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/solving_strategies/schemes/residual_based_adjoint_bossak_scheme/ResidualBasedAdjointBossakScheme/UpdateTimeSchemeAdjoints
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/solving_strategies/schemes/residual_based_adjoint_bossak_scheme/ResidualBasedAdjointBossakScheme/UpdateAuxiliaryVariable
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`static std::vector<const VariableData*> GatherVariables( const ModelPart::ElementsContainerType& rElements, std::function<void(const AdjointExtensions&, std::vector<const VariableData*>&)> GetLocalVars)`

# Calls

- [ThisThread](../../../../../functions/utilities/openmp_utils/ThisThread.md)

# Called by

- [Check](../../../../../functions/solving_strategies/schemes/residual_based_adjoint_bossak_scheme/ResidualBasedAdjointBossakScheme/Check.md)
- [UpdateTimeSchemeAdjoints](../../../../../functions/solving_strategies/schemes/residual_based_adjoint_bossak_scheme/ResidualBasedAdjointBossakScheme/UpdateTimeSchemeAdjoints.md)
- [UpdateAuxiliaryVariable](../../../../../functions/solving_strategies/schemes/residual_based_adjoint_bossak_scheme/ResidualBasedAdjointBossakScheme/UpdateAuxiliaryVariable.md)