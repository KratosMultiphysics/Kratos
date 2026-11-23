---
type: C++ Method
title: CalculateSystemContributions
resource: solving_strategies/schemes/residual_based_adjoint_bossak_scheme.h#L244-L286
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
  - target: functions/solving_strategies/schemes/residual_based_adjoint_bossak_scheme/ResidualBasedAdjointBossakScheme/CalculatePreviousTimeStepContributions
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void CalculateSystemContributions( Element& rCurrentElement, LocalSystemMatrixType& rLHS_Contribution, LocalSystemVectorType& rRHS_Contribution, Element::EquationIdVectorType& rEquationId, const ProcessInfo& rCurrentProcessInfo) override`

# Calls

- [ThisThread](../../../../../functions/utilities/openmp_utils/ThisThread.md)
- [CheckAndResizeThreadStorage](../../../../../functions/solving_strategies/schemes/residual_based_adjoint_bossak_scheme/ResidualBasedAdjointBossakScheme/CheckAndResizeThreadStorage.md)
- [CalculatePreviousTimeStepContributions](../../../../../functions/solving_strategies/schemes/residual_based_adjoint_bossak_scheme/ResidualBasedAdjointBossakScheme/CalculatePreviousTimeStepContributions.md)