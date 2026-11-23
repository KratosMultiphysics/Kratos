---
type: C++ Function
title: GlobalRetrieveGlobalPointers
resource: utilities/global_pointer_utilities.h#L280-L318
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/utilities/global_pointer_utilities/RetrieveGlobalIndexedPointers
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/mpi/tests/cpp_tests/utilities/test_pointer_communicator/KRATOS_DISTRIBUTED_TEST_CASE_IN_SUITE
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/mpi/tests/test_utilities/mpi_cpp_test_utilities/MPICppTestUtilities/SynchronizeNodes
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`static GlobalPointersVector< typename TContainerType::value_type > GlobalRetrieveGlobalPointers( const TContainerType& rContainer, const DataCommunicator& rDataCommunicator )`

# Calls

- [RetrieveGlobalIndexedPointers](../../../functions/utilities/global_pointer_utilities/RetrieveGlobalIndexedPointers.md)

# Called by

- [KRATOS_DISTRIBUTED_TEST_CASE_IN_SUITE](../../../functions/mpi/tests/cpp_tests/utilities/test_pointer_communicator/KRATOS_DISTRIBUTED_TEST_CASE_IN_SUITE.md)
- [SynchronizeNodes](../../../functions/mpi/tests/test_utilities/mpi_cpp_test_utilities/MPICppTestUtilities/SynchronizeNodes.md)