---
type: C++ Method
title: CallFunctionComponents
resource: processes/assign_scalar_field_to_entities_process.cpp#L244-L255
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

`void AssignScalarFieldToEntitiesProcess<Element>::CallFunctionComponents( const typename Element::Pointer pEntity, const double Time, double& rValue, GenericFunctionUtility& rFunction )`

# Calls

- [GetGeometry](../../../../functions/input_output/cad_json_input/CadJsonInput/GetGeometry.md)