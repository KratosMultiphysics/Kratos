---
type: C++ Method
title: FilterMPIInterfaceNodes
resource: processes/skin_detection_process.cpp#L426-L524
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/includes/data_communicator/Send
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/includes/data_communicator/Recv
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/processes/skin_detection_process/SkinDetectionProcess<TDim>/Execute
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void SkinDetectionProcess<TDim>::FilterMPIInterfaceNodes( const std::unordered_set<IndexType>& rSetNodeIdsInterface, HashMapVectorIntType& rInverseFaceMap )`

# Calls

- [Send](../../../../functions/includes/data_communicator/Send.md)
- [Recv](../../../../functions/includes/data_communicator/Recv.md)

# Called by

- [Execute](../../../../functions/processes/skin_detection_process/SkinDetectionProcess<TDim>/Execute.md)