---
type: C++ Function
title: KRATOS_DISTRIBUTED_TEST_CASE_IN_SUITE
resource: mpi/tests/cpp_tests/utilities/test_pointer_communicator.cpp#L32-L87
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/utilities/global_pointer_utilities/RetrieveGlobalIndexedPointers
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/parallel_environment/ParallelEnvironment/UnregisterDataCommunicator
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/global_pointer_utilities/LocalRetrieveGlobalPointers
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/global_pointer_utilities/GlobalRetrieveGlobalPointers
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/pointer_map_communicator/GetApplyProxy
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`KRATOS_DISTRIBUTED_TEST_CASE_IN_SUITE(PointerCommunicator, KratosMPICoreFastSuite)`

# Calls

- [RetrieveGlobalIndexedPointers](../../../../../../functions/utilities/global_pointer_utilities/RetrieveGlobalIndexedPointers.md)
- [UnregisterDataCommunicator](../../../../../../functions/sources/parallel_environment/ParallelEnvironment/UnregisterDataCommunicator.md)
- [LocalRetrieveGlobalPointers](../../../../../../functions/utilities/global_pointer_utilities/LocalRetrieveGlobalPointers.md)
- [GlobalRetrieveGlobalPointers](../../../../../../functions/utilities/global_pointer_utilities/GlobalRetrieveGlobalPointers.md)
- [GetApplyProxy](../../../../../../functions/utilities/pointer_map_communicator/GetApplyProxy.md)