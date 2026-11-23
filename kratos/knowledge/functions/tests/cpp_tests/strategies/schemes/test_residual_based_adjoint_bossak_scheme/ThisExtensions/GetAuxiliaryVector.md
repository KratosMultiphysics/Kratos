---
type: C++ Method
title: GetAuxiliaryVector
resource: tests/cpp_tests/strategies/schemes/test_residual_based_adjoint_bossak_scheme.cpp#L399-L409
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/input_output/cad_json_input/CadJsonInput/GetGeometry
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void GetAuxiliaryVector(std::size_t NodeId, std::vector<IndirectScalar<double>>& rVector, std::size_t Step) override`

# Calls

- [GetGeometry](../../../../../../../functions/input_output/cad_json_input/CadJsonInput/GetGeometry.md)