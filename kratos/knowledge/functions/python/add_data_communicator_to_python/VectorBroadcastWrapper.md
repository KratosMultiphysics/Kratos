---
type: C++ Function
title: VectorBroadcastWrapper
resource: python/add_data_communicator_to_python.cpp#L27-L58
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/includes/data_communicator/Broadcast
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`std::vector<TValue> VectorBroadcastWrapper( const DataCommunicator& rSelf, void (DataCommunicator::*pBroadcastMethod)(std::vector<TValue>&, const int) const, const std::vector<TValue>& rSourceValues, const int SourceRank)`

# Calls

- [Broadcast](../../../functions/includes/data_communicator/Broadcast.md)