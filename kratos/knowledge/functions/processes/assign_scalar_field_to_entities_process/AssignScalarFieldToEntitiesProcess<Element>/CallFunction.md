---
type: C++ Method
title: CallFunction
resource: processes/assign_scalar_field_to_entities_process.cpp#L175-L193
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/input_output/cad_json_input/CadJsonInput/GetGeometry
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/includes/node/X0
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/includes/node/Y0
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/includes/node/Z0
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void AssignScalarFieldToEntitiesProcess<Element>::CallFunction( const typename Element::Pointer pEntity, const double Time, Vector& rValue, GenericFunctionUtility& rFunction )`

# Calls

- [GetGeometry](../../../../functions/input_output/cad_json_input/CadJsonInput/GetGeometry.md)
- [X0](../../../../functions/includes/node/X0.md)
- [Y0](../../../../functions/includes/node/Y0.md)
- [Z0](../../../../functions/includes/node/Z0.md)