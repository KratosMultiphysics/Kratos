---
type: C++ Method
title: CreatePointerCommunicator
resource: processes/calculate_discontinuous_distance_to_skin_process.cpp#L979-L1010
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/utilities/global_pointer_utilities/RetrieveGlobalIndexedPointersMap
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/processes/calculate_discontinuous_distance_to_skin_process/CalculateDiscontinuousDistanceToSkinProcess<TDim>/CheckAndCorrectEdgeDistances
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`GlobalPointerCommunicator<Element>::Pointer CalculateDiscontinuousDistanceToSkinProcess<TDim>::CreatePointerCommunicator()`

# Calls

- [RetrieveGlobalIndexedPointersMap](../../../../functions/utilities/global_pointer_utilities/RetrieveGlobalIndexedPointersMap.md)

# Called by

- [CheckAndCorrectEdgeDistances](../../../../functions/processes/calculate_discontinuous_distance_to_skin_process/CalculateDiscontinuousDistanceToSkinProcess<TDim>/CheckAndCorrectEdgeDistances.md)