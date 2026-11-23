---
type: C++ Function
title: RetrieveGlobalIndexedPointers
resource: utilities/global_pointer_utilities.h#L329-L352
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/utilities/global_pointer_utilities/RetrieveGlobalIndexedPointersMap
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/mpi/tests/cpp_tests/utilities/test_pointer_communicator/KRATOS_DISTRIBUTED_TEST_CASE_IN_SUITE
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/mpi/utilities/debug_utilities/CheckNonHistoricalVariable
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/mpi/utilities/debug_utilities/CheckHistoricalVariable
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/global_pointer_utilities/LocalRetrieveGlobalPointers
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/global_pointer_utilities/GlobalRetrieveGlobalPointers
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`static GlobalPointersVector< typename TContainerType::value_type > RetrieveGlobalIndexedPointers( const TContainerType& rContainer, const std::vector<int>& rIdList, const DataCommunicator& rDataCommunicator )`

# Calls

- [RetrieveGlobalIndexedPointersMap](../../../functions/utilities/global_pointer_utilities/RetrieveGlobalIndexedPointersMap.md)

# Called by

- [KRATOS_DISTRIBUTED_TEST_CASE_IN_SUITE](../../../functions/mpi/tests/cpp_tests/utilities/test_pointer_communicator/KRATOS_DISTRIBUTED_TEST_CASE_IN_SUITE.md)
- [CheckNonHistoricalVariable](../../../functions/mpi/utilities/debug_utilities/CheckNonHistoricalVariable.md)
- [CheckHistoricalVariable](../../../functions/mpi/utilities/debug_utilities/CheckHistoricalVariable.md)
- [LocalRetrieveGlobalPointers](../../../functions/utilities/global_pointer_utilities/LocalRetrieveGlobalPointers.md)
- [GlobalRetrieveGlobalPointers](../../../functions/utilities/global_pointer_utilities/GlobalRetrieveGlobalPointers.md)