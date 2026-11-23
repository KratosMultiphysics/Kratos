---
type: C++ Method
title: GetEntityTypeName
resource: utilities/entities_utilities.cpp#L111-L122
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/utilities/entities_utilities/EntitityIdentifier<TEntity>/GenerateSingleType
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/entities_utilities/EntitityIdentifier<TEntity>/GenerateMultipleTypes
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`std::string EntitityIdentifier<TEntity>::GetEntityTypeName() const`

# Called by

- [GenerateSingleType](../../../../functions/utilities/entities_utilities/EntitityIdentifier<TEntity>/GenerateSingleType.md)
- [GenerateMultipleTypes](../../../../functions/utilities/entities_utilities/EntitityIdentifier<TEntity>/GenerateMultipleTypes.md)