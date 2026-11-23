---
type: C++ Function
title: Recv
resource: includes/data_communicator.h#L574-L577
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/includes/data_communicator/RecvImpl
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/input_output/stl_io/StlIO/WriteEntityBlockMPI
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/input_output/stl_io/StlIO/WriteGeometryBlockMPI
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/mpi/tests/cpp_tests/sources/test_mpi_data_communicator/MPIDataCommunicatorSendAndRecvIntegralTypeTest
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/mpi/tests/cpp_tests/sources/test_mpi_data_communicator/KRATOS_TEST_CASE_IN_SUITE
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/mpi/tests/cpp_tests/sources/test_mpi_serialized_communication/KRATOS_TEST_CASE_IN_SUITE
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/mpi/utilities/gather_modelpart_utility/GatherModelPartUtility/GatherEntityFromOtherPartitions
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/mpi/utilities/parallel_fill_communicator/ParallelFillCommunicator/ComputeCommunicationPlan
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/skin_detection_process/SkinDetectionProcess<TDim>/FillAuxiliaryModelPart
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/skin_detection_process/SkinDetectionProcess<TDim>/FilterMPIInterfaceNodes
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

`void Recv(TObject& rRecvObject, const int RecvSource, const int RecvTag = 0) const`

# Called by

- [RecvImpl](../../../functions/includes/data_communicator/RecvImpl.md)
- [WriteEntityBlockMPI](../../../functions/input_output/stl_io/StlIO/WriteEntityBlockMPI.md)
- [WriteGeometryBlockMPI](../../../functions/input_output/stl_io/StlIO/WriteGeometryBlockMPI.md)
- [MPIDataCommunicatorSendAndRecvIntegralTypeTest](../../../functions/mpi/tests/cpp_tests/sources/test_mpi_data_communicator/MPIDataCommunicatorSendAndRecvIntegralTypeTest.md)
- [KRATOS_TEST_CASE_IN_SUITE](../../../functions/mpi/tests/cpp_tests/sources/test_mpi_data_communicator/KRATOS_TEST_CASE_IN_SUITE.md)
- [KRATOS_TEST_CASE_IN_SUITE](../../../functions/mpi/tests/cpp_tests/sources/test_mpi_serialized_communication/KRATOS_TEST_CASE_IN_SUITE.md)
- [GatherEntityFromOtherPartitions](../../../functions/mpi/utilities/gather_modelpart_utility/GatherModelPartUtility/GatherEntityFromOtherPartitions.md)
- [ComputeCommunicationPlan](../../../functions/mpi/utilities/parallel_fill_communicator/ParallelFillCommunicator/ComputeCommunicationPlan.md)
- [FillAuxiliaryModelPart](../../../functions/processes/skin_detection_process/SkinDetectionProcess<TDim>/FillAuxiliaryModelPart.md)
- [FilterMPIInterfaceNodes](../../../functions/processes/skin_detection_process/SkinDetectionProcess<TDim>/FilterMPIInterfaceNodes.md)
- [SetParallelModelPartAndSubModelPartCollectionsAndCombinations](../../../functions/utilities/assign_unique_model_part_collection_tag_utility/AssignUniqueModelPartCollectionTagUtility/SetParallelModelPartAndSubModelPartCollectionsAndCombinations.md)
- [RetrieveGlobalIndexedPointersMap](../../../functions/utilities/global_pointer_utilities/RetrieveGlobalIndexedPointersMap.md)