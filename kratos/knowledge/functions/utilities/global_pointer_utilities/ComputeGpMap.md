---
type: C++ Function
title: ComputeGpMap
resource: utilities/global_pointer_utilities.h#L486-L519
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/utilities/global_pointer_utilities/RetrieveGlobalIndexedPointersMap
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`static std::unordered_map< int, GlobalPointer<typename TContainerType::value_type> > ComputeGpMap( const TContainerType& rContainer, const std::vector<int>& rIds, const DataCommunicator& rDataCommunicator)`

# Called by

- [RetrieveGlobalIndexedPointersMap](../../../functions/utilities/global_pointer_utilities/RetrieveGlobalIndexedPointersMap.md)