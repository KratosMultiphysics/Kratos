---
type: C++ Method
title: StreamSerializer
resource: includes/stream_serializer.h#L39-L41
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/includes/serializer/pGetBuffer
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`explicit StreamSerializer(TraceType const& rTrace=SERIALIZER_NO_TRACE) : Serializer(new std::stringstream(std::ios::binary|std::ios::in|std::ios::out),rTrace)`

# Calls

- [pGetBuffer](../../../../functions/includes/serializer/pGetBuffer.md)