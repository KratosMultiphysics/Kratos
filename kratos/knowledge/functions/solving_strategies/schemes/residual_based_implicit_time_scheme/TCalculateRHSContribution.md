---
type: C++ Function
title: TCalculateRHSContribution
resource: solving_strategies/schemes/residual_based_implicit_time_scheme.h#L490-L512
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/utilities/openmp_utils/ThisThread
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/solving_strategies/schemes/residual_based_implicit_time_scheme/ResidualBasedImplicitTimeScheme/CalculateRHSContribution
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void TCalculateRHSContribution( TObjectType& rObject, LocalSystemVectorType& rRHSContribution, Element::EquationIdVectorType& rEquationId, const ProcessInfo& rCurrentProcessInfo )`

# Calls

- [ThisThread](../../../../functions/utilities/openmp_utils/ThisThread.md)

# Called by

- [CalculateRHSContribution](../../../../functions/solving_strategies/schemes/residual_based_implicit_time_scheme/ResidualBasedImplicitTimeScheme/CalculateRHSContribution.md)