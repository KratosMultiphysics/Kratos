---
type: C++ Function
title: AddSerializerToPython
resource: python/add_serializer_to_python.cpp#L59-L143
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/includes/stream_serializer/StreamSerializer/GetStringRepresentation
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/includes/serializer/GetTraceType
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/python/kratos_python/PYBIND11_MODULE
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void AddSerializerToPython(pybind11::module& m)`

# Calls

- [GetStringRepresentation](../../../functions/includes/stream_serializer/StreamSerializer/GetStringRepresentation.md)
- [GetTraceType](../../../functions/includes/serializer/GetTraceType.md)

# Called by

- [PYBIND11_MODULE](../../../functions/python/kratos_python/PYBIND11_MODULE.md)