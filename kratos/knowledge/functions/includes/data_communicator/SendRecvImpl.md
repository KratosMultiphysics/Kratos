---
type: C++ Function
title: SendRecvImpl
resource: includes/data_communicator.h#L860-L866
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/includes/stream_serializer/StreamSerializer/GetStringRepresentation
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`virtual void SendRecvImpl( const std::string& rSendValues, const int SendDestination, const int SendTag, std::string& rRecvValues, const int RecvSource, const int RecvTag) const`

# Calls

- [GetStringRepresentation](../../../functions/includes/stream_serializer/StreamSerializer/GetStringRepresentation.md)