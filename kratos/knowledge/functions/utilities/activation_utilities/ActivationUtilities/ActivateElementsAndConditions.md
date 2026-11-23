---
type: C++ Method
title: ActivateElementsAndConditions
resource: utilities/activation_utilities.h#L103-L169
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/input_output/cad_json_input/CadJsonInput/GetGeometry
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/containers/flags/Flip
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void ActivateElementsAndConditions( ModelPart& rmodel_part, const Variable< double >& rVariable, const double reference_value, bool active_if_lower_than_reference)`

# Calls

- [GetGeometry](../../../../functions/input_output/cad_json_input/CadJsonInput/GetGeometry.md)
- [Flip](../../../../functions/containers/flags/Flip.md)