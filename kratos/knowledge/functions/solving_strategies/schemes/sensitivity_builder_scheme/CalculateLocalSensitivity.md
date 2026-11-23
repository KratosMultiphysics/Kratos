---
type: C++ Function
title: CalculateLocalSensitivity
resource: solving_strategies/schemes/sensitivity_builder_scheme.h#L587-L621
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
  - target: functions/solving_strategies/schemes/sensitivity_builder_scheme/CalculateLocalSensitivityAndGlobalPointersVector
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void CalculateLocalSensitivity( TEntityType& rCurrentEntity, AdjointResponseFunction& rResponseFunction, Vector& rSensitivity, const Variable<TDataType>& rVariable, const ProcessInfo& rCurrentProcessInfo)`

# Calls

- [ThisThread](../../../../functions/utilities/openmp_utils/ThisThread.md)
- [prod](../../../../functions/includes/amatrix_interface/prod.md)

# Called by

- [CalculateLocalSensitivityAndGlobalPointersVector](../../../../functions/solving_strategies/schemes/sensitivity_builder_scheme/CalculateLocalSensitivityAndGlobalPointersVector.md)