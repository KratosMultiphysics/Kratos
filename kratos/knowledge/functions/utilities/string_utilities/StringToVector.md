---
type: C++ Function
title: StringToVector
resource: utilities/string_utilities.cpp#L263-L308
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/utilities/string_utilities/Trim
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/containers/pointer_vector_set/PointerVectorSet/pop_back
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`std::vector<TType> StringToVector(const std::string& rInputString)`

# Calls

- [Trim](../../../functions/utilities/string_utilities/Trim.md)
- [pop_back](../../../functions/containers/pointer_vector_set/PointerVectorSet/pop_back.md)