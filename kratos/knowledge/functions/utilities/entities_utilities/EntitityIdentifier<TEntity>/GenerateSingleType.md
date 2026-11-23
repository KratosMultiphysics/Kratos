---
type: C++ Method
title: GenerateSingleType
resource: utilities/entities_utilities.cpp#L128-L140
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/utilities/entities_utilities/EntitityIdentifier<TEntity>/GetEntityTypeName
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/input_output/cad_json_input/CadJsonInput/GetGeometry
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/utilities/entities_utilities/EntitityIdentifier<TEntity>/EntitityIdentifier
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void EntitityIdentifier<TEntity>::GenerateSingleType(const std::string& rName)`

# Calls

- [GetEntityTypeName](../../../../functions/utilities/entities_utilities/EntitityIdentifier<TEntity>/GetEntityTypeName.md)
- [GetGeometry](../../../../functions/input_output/cad_json_input/CadJsonInput/GetGeometry.md)

# Called by

- [EntitityIdentifier](../../../../functions/utilities/entities_utilities/EntitityIdentifier<TEntity>/EntitityIdentifier.md)