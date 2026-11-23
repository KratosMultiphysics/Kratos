---
type: C++ Function
title: AddInterfaceToAccessorFold
resource: python/add_properties_to_python.cpp#L152-L171
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/includes/properties/SetAccessor
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/includes/properties/HasAccessor
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void AddInterfaceToAccessorFold(pybind11::class_<Properties, Properties::Pointer, Properties::BaseType>& pyProperties)`

# Calls

- [SetAccessor](../../../functions/includes/properties/SetAccessor.md)
- [HasAccessor](../../../functions/includes/properties/HasAccessor.md)