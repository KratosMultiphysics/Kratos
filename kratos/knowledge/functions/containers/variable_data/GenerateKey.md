---
type: C++ Function
title: GenerateKey
resource: containers/variable_data.h#L285-L301
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/containers/variable_data/VariableData
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`static constexpr KeyType GenerateKey( const std::string_view& Name, const std::size_t Size, const bool IsComponent, const char ComponentIndex)`

# Called by

- [VariableData](../../../functions/containers/variable_data/VariableData.md)