---
type: C++ Method
title: GenerateTemplatedTypes
resource: utilities/entities_utilities.cpp#L168-L184
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/utilities/string_utilities/ReplaceAllSubstrings
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

`void EntitityIdentifier<TEntity>::GenerateTemplatedTypes(const std::string& rName)`

# Calls

- [ReplaceAllSubstrings](../../../../functions/utilities/string_utilities/ReplaceAllSubstrings.md)
- [GetGeometry](../../../../functions/input_output/cad_json_input/CadJsonInput/GetGeometry.md)

# Called by

- [EntitityIdentifier](../../../../functions/utilities/entities_utilities/EntitityIdentifier<TEntity>/EntitityIdentifier.md)