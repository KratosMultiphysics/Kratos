---
type: C++ Function
title: AddRegistryToPython
resource: python/add_registry_to_python.cpp#L52-L78
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/sources/registry_item/RegistryItem/KeyConstBegin
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/registry_item/RegistryItem/KeyConstEnd
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/python/kratos_python/PYBIND11_MODULE
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void AddRegistryToPython(pybind11::module& m)`

# Calls

- [KeyConstBegin](../../../functions/sources/registry_item/RegistryItem/KeyConstBegin.md)
- [KeyConstEnd](../../../functions/sources/registry_item/RegistryItem/KeyConstEnd.md)

# Called by

- [PYBIND11_MODULE](../../../functions/python/kratos_python/PYBIND11_MODULE.md)