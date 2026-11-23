---
type: C++ Function
title: CheckNonHistoricalVariable
resource: mpi/utilities/debug_utilities.h#L95-L114
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/utilities/global_pointer_utilities/RetrieveGlobalIndexedPointersMap
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/global_pointer_utilities/RetrieveGlobalIndexedPointers
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`static void CheckNonHistoricalVariable( ModelPart & rModelPart, const TContainerType & rContainer, const TVarType & rVariable)`

# Calls

- [RetrieveGlobalIndexedPointersMap](../../../../functions/utilities/global_pointer_utilities/RetrieveGlobalIndexedPointersMap.md)
- [RetrieveGlobalIndexedPointers](../../../../functions/utilities/global_pointer_utilities/RetrieveGlobalIndexedPointers.md)