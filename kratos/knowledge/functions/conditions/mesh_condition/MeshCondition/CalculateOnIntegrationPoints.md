---
type: C++ Method
title: CalculateOnIntegrationPoints
resource: conditions/mesh_condition.cpp#L170-L182
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/input_output/cad_json_input/CadJsonInput/GetGeometry
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/geometries/geometry/GetDefaultIntegrationMethod
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void MeshCondition::CalculateOnIntegrationPoints( const Variable<bool>& rVariable, std::vector<bool>& rOutput, const ProcessInfo& rCurrentProcessInfo )`

# Calls

- [GetGeometry](../../../../functions/input_output/cad_json_input/CadJsonInput/GetGeometry.md)
- [GetDefaultIntegrationMethod](../../../../functions/geometries/geometry/GetDefaultIntegrationMethod.md)