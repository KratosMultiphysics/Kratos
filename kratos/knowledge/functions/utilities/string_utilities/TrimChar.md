---
type: C++ Function
title: TrimChar
resource: utilities/string_utilities.cpp#L192-L203
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/utilities/string_utilities/TrimLeft
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/string_utilities/TrimRight
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`std::function<bool(std::string::value_type)> TrimChar(const bool RemoveNullChar)`

# Called by

- [TrimLeft](../../../functions/utilities/string_utilities/TrimLeft.md)
- [TrimRight](../../../functions/utilities/string_utilities/TrimRight.md)