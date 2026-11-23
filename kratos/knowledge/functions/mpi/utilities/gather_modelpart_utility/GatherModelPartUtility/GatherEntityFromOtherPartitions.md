---
type: C++ Method
title: GatherEntityFromOtherPartitions
resource: mpi/utilities/gather_modelpart_utility.cpp#L301-L461
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/includes/data_communicator/Recv
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/includes/data_communicator/Send
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/includes/serializer/pGetBuffer
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void GatherModelPartUtility::GatherEntityFromOtherPartitions( ModelPart& rModelPart, const std::map<int, std::vector<std::size_t>>& rEntitiesToBring, const int EchoLevel )`

# Calls

- [Recv](../../../../../functions/includes/data_communicator/Recv.md)
- [Send](../../../../../functions/includes/data_communicator/Send.md)
- [pGetBuffer](../../../../../functions/includes/serializer/pGetBuffer.md)