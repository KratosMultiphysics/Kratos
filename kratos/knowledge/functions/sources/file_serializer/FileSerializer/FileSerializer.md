---
type: C++ Method
title: FileSerializer
resource: sources/file_serializer.cpp#L22-L33
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/includes/serializer/SetBuffer
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/includes/serializer/pGetBuffer
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`FileSerializer::FileSerializer(std::string const& Filename, Serializer::TraceType const& rTrace, const bool DataOnly) : Serializer(nullptr, rTrace, DataOnly)`

# Calls

- [SetBuffer](../../../../functions/includes/serializer/SetBuffer.md)
- [pGetBuffer](../../../../functions/includes/serializer/pGetBuffer.md)