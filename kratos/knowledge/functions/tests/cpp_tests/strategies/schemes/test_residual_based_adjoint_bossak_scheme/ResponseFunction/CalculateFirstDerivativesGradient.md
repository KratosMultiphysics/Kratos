---
type: C++ Method
title: CalculateFirstDerivativesGradient
resource: tests/cpp_tests/strategies/schemes/test_residual_based_adjoint_bossak_scheme.cpp#L751-L759
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/input_output/cad_json_input/CadJsonInput/GetGeometry
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void CalculateFirstDerivativesGradient(const Element& rAdjointElement, const Matrix& rResidualGradient, Vector& rResponseGradient, const ProcessInfo& rProcessInfo) override`

# Calls

- [GetGeometry](../../../../../../../functions/input_output/cad_json_input/CadJsonInput/GetGeometry.md)