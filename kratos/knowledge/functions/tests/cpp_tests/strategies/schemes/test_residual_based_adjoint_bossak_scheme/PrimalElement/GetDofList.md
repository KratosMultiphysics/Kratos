---
type: C++ Method
title: GetDofList
resource: tests/cpp_tests/strategies/schemes/test_residual_based_adjoint_bossak_scheme.cpp#L200-L205
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/input_output/cad_json_input/CadJsonInput/GetGeometry
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void GetDofList(DofsVectorType& rElementalDofList, const ProcessInfo& rCurrentProcessInfo) const override`

# Calls

- [GetGeometry](../../../../../../../functions/input_output/cad_json_input/CadJsonInput/GetGeometry.md)