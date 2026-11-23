---
type: C++ Method
title: GlobalNumberOfNodes
resource: sources/communicator.cpp#L101-L104
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/input_output/vtk_output/VtkOutput/PrintOutput
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/mpi/tests/cpp_tests/sources/test_mpi_communicator/KRATOS_TEST_CASE_IN_SUITE
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/mpi/tests/cpp_tests/utilities/test_distributed_model_part_initializer/KRATOS_TEST_CASE_IN_SUITE
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/mpi/tests/test_distributed_import_model_part_utility/TestDistributedImportModelPartUtility/__execute_test
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/mpi/tests/test_distributed_model_part_initializer/TestDistributedModelPartInitializer/test_with_real_model
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/mpi/tests/test_mpi_communicator/TestMPICommunicator/test_GlobalNumberOf_Methods
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/calculate_discontinuous_distance_to_skin_process/CalculateDiscontinuousDistanceToSkinProcess<TDim>/CalculateCharacteristicLength
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/edge_based_gradient_recovery_process/EdgeBasedGradientRecoveryProcess/CheckAndAssignSettings
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_combine_model_part_modeler/TestCombineModelPartModeler/test_CombineModelParts
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`Communicator::SizeType Communicator::GlobalNumberOfNodes() const`

# Called by

- [PrintOutput](../../../../functions/input_output/vtk_output/VtkOutput/PrintOutput.md)
- [KRATOS_TEST_CASE_IN_SUITE](../../../../functions/mpi/tests/cpp_tests/sources/test_mpi_communicator/KRATOS_TEST_CASE_IN_SUITE.md)
- [KRATOS_TEST_CASE_IN_SUITE](../../../../functions/mpi/tests/cpp_tests/utilities/test_distributed_model_part_initializer/KRATOS_TEST_CASE_IN_SUITE.md)
- [__execute_test](../../../../functions/mpi/tests/test_distributed_import_model_part_utility/TestDistributedImportModelPartUtility/__execute_test.md)
- [test_with_real_model](../../../../functions/mpi/tests/test_distributed_model_part_initializer/TestDistributedModelPartInitializer/test_with_real_model.md)
- [test_GlobalNumberOf_Methods](../../../../functions/mpi/tests/test_mpi_communicator/TestMPICommunicator/test_GlobalNumberOf_Methods.md)
- [CalculateCharacteristicLength](../../../../functions/processes/calculate_discontinuous_distance_to_skin_process/CalculateDiscontinuousDistanceToSkinProcess<TDim>/CalculateCharacteristicLength.md)
- [CheckAndAssignSettings](../../../../functions/processes/edge_based_gradient_recovery_process/EdgeBasedGradientRecoveryProcess/CheckAndAssignSettings.md)
- [test_CombineModelParts](../../../../functions/tests/test_combine_model_part_modeler/TestCombineModelPartModeler/test_CombineModelParts.md)