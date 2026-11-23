---
type: C++ Function
title: SetValuesOnIntegrationPointsVector
resource: python/add_mesh_to_python.cpp#L212-L225
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/input_output/cad_json_input/CadJsonInput/GetGeometry
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void SetValuesOnIntegrationPointsVector( TObject& dummy, const Variable<Vector>& rVariable, pybind11::list values_list, unsigned int len_values_list_item, const ProcessInfo& rCurrentProcessInfo )`

# Calls

- [GetGeometry](../../../functions/input_output/cad_json_input/CadJsonInput/GetGeometry.md)