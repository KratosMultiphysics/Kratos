---
type: C++ Method
title: GetDofList
resource: sources/periodic_condition.cpp#L100-L119
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/input_output/cad_json_input/CadJsonInput/GetGeometry
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/geometries/geometry/PointsNumber
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/containers/periodic_variables_container/PeriodicVariablesContainer/DoubleVariablesBegin
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/containers/periodic_variables_container/PeriodicVariablesContainer/DoubleVariablesEnd
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void PeriodicCondition::GetDofList(DofsVectorType& ElementalDofList, const ProcessInfo& rCurrentProcessInfo) const`

# Calls

- [GetGeometry](../../../../functions/input_output/cad_json_input/CadJsonInput/GetGeometry.md)
- [PointsNumber](../../../../functions/geometries/geometry/PointsNumber.md)
- [DoubleVariablesBegin](../../../../functions/containers/periodic_variables_container/PeriodicVariablesContainer/DoubleVariablesBegin.md)
- [DoubleVariablesEnd](../../../../functions/containers/periodic_variables_container/PeriodicVariablesContainer/DoubleVariablesEnd.md)