---
type: C++ Function
title: BroadcastImpl
resource: includes/data_communicator.h#L817-L849
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/includes/stream_serializer/StreamSerializer/GetStringRepresentation
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/includes/data_communicator/Broadcast
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void BroadcastImpl(TObject& rBroadcastObject, const int SourceRank) const`

# Calls

- [GetStringRepresentation](../../../functions/includes/stream_serializer/StreamSerializer/GetStringRepresentation.md)
- [Broadcast](../../../functions/includes/data_communicator/Broadcast.md)