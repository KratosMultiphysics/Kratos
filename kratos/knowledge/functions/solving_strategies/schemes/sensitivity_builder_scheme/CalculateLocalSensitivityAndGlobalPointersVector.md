---
type: C++ Function
title: CalculateLocalSensitivityAndGlobalPointersVector
resource: solving_strategies/schemes/sensitivity_builder_scheme.h#L537-L557
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/solving_strategies/schemes/sensitivity_builder_scheme/CalculateLocalSensitivity
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/input_output/cad_json_input/CadJsonInput/GetGeometry
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/geometries/geometry/PointsNumber
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void CalculateLocalSensitivityAndGlobalPointersVector( TEntityType& rEntity, AdjointResponseFunction& rResponseFunction, Vector& rSensitivityVector, GlobalPointersVector<TDerivativeEntityType>& rGPSensitivityVector, const Variable<TDataType>& rVariable, const ProcessInfo& rProcessInfo)`

# Calls

- [CalculateLocalSensitivity](../../../../functions/solving_strategies/schemes/sensitivity_builder_scheme/CalculateLocalSensitivity.md)
- [GetGeometry](../../../../functions/input_output/cad_json_input/CadJsonInput/GetGeometry.md)
- [PointsNumber](../../../../functions/geometries/geometry/PointsNumber.md)