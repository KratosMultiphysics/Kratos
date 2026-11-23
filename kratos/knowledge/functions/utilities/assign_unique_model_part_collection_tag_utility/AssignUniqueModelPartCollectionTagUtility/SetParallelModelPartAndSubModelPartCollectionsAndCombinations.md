---
type: C++ Method
title: SetParallelModelPartAndSubModelPartCollectionsAndCombinations
resource: utilities/assign_unique_model_part_collection_tag_utility.cpp#L175-L231
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
  - target: functions/includes/data_communicator/Broadcast
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/utilities/assign_unique_model_part_collection_tag_utility/AssignUniqueModelPartCollectionTagUtility/ComputeTags
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void AssignUniqueModelPartCollectionTagUtility::SetParallelModelPartAndSubModelPartCollectionsAndCombinations( IndexStringMapType& rCollections, IndexSetIndexMapType& rCombinations, IndexType& rTag)`

# Calls

- [Send](../../../../functions/includes/data_communicator/Send.md)
- [Recv](../../../../functions/includes/data_communicator/Recv.md)
- [Broadcast](../../../../functions/includes/data_communicator/Broadcast.md)

# Called by

- [ComputeTags](../../../../functions/utilities/assign_unique_model_part_collection_tag_utility/AssignUniqueModelPartCollectionTagUtility/ComputeTags.md)