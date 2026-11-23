---
type: C++ Method
title: AssignTimeDependentValue
resource: processes/assign_scalar_field_to_entities_process.cpp#L444-L461
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/input_output/cad_json_input/CadJsonInput/GetGeometry
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void AssignScalarFieldToEntitiesProcess<Condition>::AssignTimeDependentValue( const typename Condition::Pointer pEntity, const double Time, Vector& rValue, const double Value )`

# Calls

- [GetGeometry](../../../../functions/input_output/cad_json_input/CadJsonInput/GetGeometry.md)