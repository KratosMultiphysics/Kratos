---
type: C++ Function
title: SetValuesOnIntegrationPointsConstitutiveLaw
resource: python/add_mesh_to_python.cpp#L235-L247
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/input_output/cad_json_input/CadJsonInput/GetGeometry
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void SetValuesOnIntegrationPointsConstitutiveLaw( Element& dummy, const Variable<ConstitutiveLaw::Pointer>& rVariable, pybind11::list values_list, const ProcessInfo& rCurrentProcessInfo )`

# Calls

- [GetGeometry](../../../functions/input_output/cad_json_input/CadJsonInput/GetGeometry.md)