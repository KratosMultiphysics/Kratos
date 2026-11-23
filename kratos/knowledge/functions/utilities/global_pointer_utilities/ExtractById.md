---
type: C++ Function
title: ExtractById
resource: utilities/global_pointer_utilities.h#L466-L476
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/utilities/global_pointer_utilities/RetrieveGlobalIndexedPointersMap
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`static std::unordered_map< int, GPType > ExtractById( std::unordered_map< int, GPType >& rGPList, const std::vector<int>& rIds)`

# Called by

- [RetrieveGlobalIndexedPointersMap](../../../functions/utilities/global_pointer_utilities/RetrieveGlobalIndexedPointersMap.md)