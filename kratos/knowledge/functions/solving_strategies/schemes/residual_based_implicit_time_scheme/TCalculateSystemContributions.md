---
type: C++ Function
title: TCalculateSystemContributions
resource: solving_strategies/schemes/residual_based_implicit_time_scheme.h#L455-L480
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/utilities/openmp_utils/ThisThread
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/solving_strategies/schemes/residual_based_implicit_time_scheme/ResidualBasedImplicitTimeScheme/CalculateSystemContributions
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void TCalculateSystemContributions( TObjectType& rObject, LocalSystemMatrixType& rLHSContribution, LocalSystemVectorType& rRHSContribution, Element::EquationIdVectorType& EquationId, const ProcessInfo& rCurrentProcessInfo )`

# Calls

- [ThisThread](../../../../functions/utilities/openmp_utils/ThisThread.md)

# Called by

- [CalculateSystemContributions](../../../../functions/solving_strategies/schemes/residual_based_implicit_time_scheme/ResidualBasedImplicitTimeScheme/CalculateSystemContributions.md)