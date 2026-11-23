---
type: C++ Function
title: Broadcast
resource: includes/data_communicator.h#L474-L477
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/includes/data_communicator/BroadcastImpl
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/mpi/sources/mpi_data_communicator/MPIDataCommunicator/ScatterDetail
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/mpi/tests/cpp_tests/sources/test_data_communicator/KRATOS_DISTRIBUTED_TEST_CASE_IN_SUITE
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/mpi/tests/cpp_tests/sources/test_mpi_data_communicator/MPIDataCommunicatorBroadcastIntegralTypeTest
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/mpi/tests/cpp_tests/sources/test_mpi_data_communicator/KRATOS_TEST_CASE_IN_SUITE
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/mpi/tests/cpp_tests/sources/test_mpi_data_communicator/MPIDataCommunicatorBroadcastIntegralTypeVectorTest
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/mpi/tests/cpp_tests/sources/test_mpi_serialized_communication/KRATOS_TEST_CASE_IN_SUITE
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/mpi/tests/test_mpi_data_communicator_python/TestMPIDataCommunicatorPython/testBroadcastOperations
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/mpi/utilities/distributed_model_part_initializer/DistributedModelPartInitializer/CopySubModelPartStructure
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/mpi/utilities/parallel_fill_communicator/ParallelFillCommunicator/ComputeCommunicationPlan
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python/add_data_communicator_to_python/VectorBroadcastWrapper
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python/add_data_communicator_to_python/AddDataCommunicatorMethodForDataType
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python/add_data_communicator_to_python/AddDataCommunicatorToPython
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/assign_unique_model_part_collection_tag_utility/AssignUniqueModelPartCollectionTagUtility/SetParallelModelPartAndSubModelPartCollectionsAndCombinations
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/global_pointer_utilities/RetrieveGlobalIndexedPointersMap
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void Broadcast(TObject& rBroadcastObject, const int SourceRank) const`

# Called by

- [BroadcastImpl](../../../functions/includes/data_communicator/BroadcastImpl.md)
- [ScatterDetail](../../../functions/mpi/sources/mpi_data_communicator/MPIDataCommunicator/ScatterDetail.md)
- [KRATOS_DISTRIBUTED_TEST_CASE_IN_SUITE](../../../functions/mpi/tests/cpp_tests/sources/test_data_communicator/KRATOS_DISTRIBUTED_TEST_CASE_IN_SUITE.md)
- [MPIDataCommunicatorBroadcastIntegralTypeTest](../../../functions/mpi/tests/cpp_tests/sources/test_mpi_data_communicator/MPIDataCommunicatorBroadcastIntegralTypeTest.md)
- [KRATOS_TEST_CASE_IN_SUITE](../../../functions/mpi/tests/cpp_tests/sources/test_mpi_data_communicator/KRATOS_TEST_CASE_IN_SUITE.md)
- [MPIDataCommunicatorBroadcastIntegralTypeVectorTest](../../../functions/mpi/tests/cpp_tests/sources/test_mpi_data_communicator/MPIDataCommunicatorBroadcastIntegralTypeVectorTest.md)
- [KRATOS_TEST_CASE_IN_SUITE](../../../functions/mpi/tests/cpp_tests/sources/test_mpi_serialized_communication/KRATOS_TEST_CASE_IN_SUITE.md)
- [testBroadcastOperations](../../../functions/mpi/tests/test_mpi_data_communicator_python/TestMPIDataCommunicatorPython/testBroadcastOperations.md)
- [CopySubModelPartStructure](../../../functions/mpi/utilities/distributed_model_part_initializer/DistributedModelPartInitializer/CopySubModelPartStructure.md)
- [ComputeCommunicationPlan](../../../functions/mpi/utilities/parallel_fill_communicator/ParallelFillCommunicator/ComputeCommunicationPlan.md)
- [VectorBroadcastWrapper](../../../functions/python/add_data_communicator_to_python/VectorBroadcastWrapper.md)
- [AddDataCommunicatorMethodForDataType](../../../functions/python/add_data_communicator_to_python/AddDataCommunicatorMethodForDataType.md)
- [AddDataCommunicatorToPython](../../../functions/python/add_data_communicator_to_python/AddDataCommunicatorToPython.md)
- [SetParallelModelPartAndSubModelPartCollectionsAndCombinations](../../../functions/utilities/assign_unique_model_part_collection_tag_utility/AssignUniqueModelPartCollectionTagUtility/SetParallelModelPartAndSubModelPartCollectionsAndCombinations.md)
- [RetrieveGlobalIndexedPointersMap](../../../functions/utilities/global_pointer_utilities/RetrieveGlobalIndexedPointersMap.md)