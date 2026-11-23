---
type: C++ Function
title: InvertNormalForFlag
resource: utilities/mortar_utilities.h#L252-L271
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/containers/flags/IsDefined
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/input_output/cad_json_input/CadJsonInput/GetGeometry
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/utilities/mortar_utilities/InvertNormal
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void InvertNormalForFlag( TContainerType& rContainer, const Flags Flag )`

# Calls

- [IsDefined](../../../functions/containers/flags/IsDefined.md)
- [GetGeometry](../../../functions/input_output/cad_json_input/CadJsonInput/GetGeometry.md)

# Called by

- [InvertNormal](../../../functions/utilities/mortar_utilities/InvertNormal.md)