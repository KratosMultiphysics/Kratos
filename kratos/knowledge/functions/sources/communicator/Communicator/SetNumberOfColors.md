---
type: C++ Method
title: SetNumberOfColors
resource: sources/communicator.cpp#L126-L144
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/modeler/combine_model_part_modeler/CombineModelPartModeler/CreateCommunicators
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/modeler/combine_model_part_modeler/CombineModelPartModeler/DuplicateCommunicatorData
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/modeler/connectivity_preserve_modeler/ConnectivityPreserveModeler/DuplicateCommunicatorData
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/mpi/utilities/parallel_fill_communicator/ParallelFillCommunicator/InitializeParallelCommunicationMeshes
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/mpi/utilities/parallel_fill_communicator/ParallelFillCommunicator/ComputeCommunicationPlanSubModelPart
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/model_part_io/ModelPartIO/ReadCommunicatorDataBlock
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/sources/test_communicator/KRATOS_TEST_CASE_IN_SUITE
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/sources/test_model_part_io/KRATOS_TEST_CASE_IN_SUITE
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/model_part_operation_utilities/ModelPartOperationUtilities/SetCommunicator
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void Communicator::SetNumberOfColors(SizeType NewNumberOfColors)`

# Called by

- [CreateCommunicators](../../../../functions/modeler/combine_model_part_modeler/CombineModelPartModeler/CreateCommunicators.md)
- [DuplicateCommunicatorData](../../../../functions/modeler/combine_model_part_modeler/CombineModelPartModeler/DuplicateCommunicatorData.md)
- [DuplicateCommunicatorData](../../../../functions/modeler/connectivity_preserve_modeler/ConnectivityPreserveModeler/DuplicateCommunicatorData.md)
- [InitializeParallelCommunicationMeshes](../../../../functions/mpi/utilities/parallel_fill_communicator/ParallelFillCommunicator/InitializeParallelCommunicationMeshes.md)
- [ComputeCommunicationPlanSubModelPart](../../../../functions/mpi/utilities/parallel_fill_communicator/ParallelFillCommunicator/ComputeCommunicationPlanSubModelPart.md)
- [ReadCommunicatorDataBlock](../../../../functions/sources/model_part_io/ModelPartIO/ReadCommunicatorDataBlock.md)
- [KRATOS_TEST_CASE_IN_SUITE](../../../../functions/tests/cpp_tests/sources/test_communicator/KRATOS_TEST_CASE_IN_SUITE.md)
- [KRATOS_TEST_CASE_IN_SUITE](../../../../functions/tests/cpp_tests/sources/test_model_part_io/KRATOS_TEST_CASE_IN_SUITE.md)
- [SetCommunicator](../../../../functions/utilities/model_part_operation_utilities/ModelPartOperationUtilities/SetCommunicator.md)