---
type: C++ Method
title: CalculateSensitivityMatrix
resource: tests/cpp_tests/strategies/schemes/test_residual_based_adjoint_bossak_scheme.cpp#L518-L537
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/input_output/cad_json_input/CadJsonInput/GetGeometry
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void CalculateSensitivityMatrix(const Variable<double>& rDesignVariable, Matrix& rOutput, const ProcessInfo& rCurrentProcessInfo) override`

# Calls

- [GetGeometry](../../../../../../../functions/input_output/cad_json_input/CadJsonInput/GetGeometry.md)