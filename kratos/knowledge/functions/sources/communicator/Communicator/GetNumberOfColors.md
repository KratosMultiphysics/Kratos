---
type: C++ Method
title: GetNumberOfColors
resource: sources/communicator.cpp#L121-L124
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
  - target: functions/mpi/tests/cpp_tests/sources/test_mpi_communicator/KRATOS_TEST_CASE_IN_SUITE
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/mpi/tests/test_mpi_model_part/TestMPIModelPart/test_remove_nodes_parallel_interfaces
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/mpi/utilities/mpi_normal_calculation_utilities/MPINormalCalculationUtils/InitializeNormalData
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/model_part_io/ModelPartIO/ReadCommunicatorLocalNodesBlock
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/model_part_io/ModelPartIO/ReadCommunicatorGhostNodesBlock
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/sources/test_communicator/KRATOS_TEST_CASE_IN_SUITE
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/model_part_operation_utilities/ModelPartOperationUtilities/SetCommunicator
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`Communicator::SizeType Communicator::GetNumberOfColors() const`

# Called by

- [CreateCommunicators](../../../../functions/modeler/combine_model_part_modeler/CombineModelPartModeler/CreateCommunicators.md)
- [DuplicateCommunicatorData](../../../../functions/modeler/combine_model_part_modeler/CombineModelPartModeler/DuplicateCommunicatorData.md)
- [DuplicateCommunicatorData](../../../../functions/modeler/connectivity_preserve_modeler/ConnectivityPreserveModeler/DuplicateCommunicatorData.md)
- [KRATOS_TEST_CASE_IN_SUITE](../../../../functions/mpi/tests/cpp_tests/sources/test_mpi_communicator/KRATOS_TEST_CASE_IN_SUITE.md)
- [test_remove_nodes_parallel_interfaces](../../../../functions/mpi/tests/test_mpi_model_part/TestMPIModelPart/test_remove_nodes_parallel_interfaces.md)
- [InitializeNormalData](../../../../functions/mpi/utilities/mpi_normal_calculation_utilities/MPINormalCalculationUtils/InitializeNormalData.md)
- [ReadCommunicatorLocalNodesBlock](../../../../functions/sources/model_part_io/ModelPartIO/ReadCommunicatorLocalNodesBlock.md)
- [ReadCommunicatorGhostNodesBlock](../../../../functions/sources/model_part_io/ModelPartIO/ReadCommunicatorGhostNodesBlock.md)
- [KRATOS_TEST_CASE_IN_SUITE](../../../../functions/tests/cpp_tests/sources/test_communicator/KRATOS_TEST_CASE_IN_SUITE.md)
- [SetCommunicator](../../../../functions/utilities/model_part_operation_utilities/ModelPartOperationUtilities/SetCommunicator.md)