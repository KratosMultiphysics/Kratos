---
type: C++ Method
title: GlobalNumberOfElements
resource: sources/communicator.cpp#L106-L109
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/input_output/ensight_output/EnSightOutput/EnSightOutput
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/input_output/ensight_output/EnSightOutput/GetEntityType
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/input_output/vtk_output/VtkOutput/VtkOutput
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/input_output/vtk_output/VtkOutput/GetEntityType
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/input_output/vtk_output/VtkOutput/PrintOutput
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/mpi/tests/cpp_tests/sources/test_mpi_communicator/KRATOS_TEST_CASE_IN_SUITE
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
  - target: functions/tests/test_combine_model_part_modeler/TestCombineModelPartModeler/test_CombineModelParts
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`Communicator::SizeType Communicator::GlobalNumberOfElements() const`

# Called by

- [EnSightOutput](../../../../functions/input_output/ensight_output/EnSightOutput/EnSightOutput.md)
- [GetEntityType](../../../../functions/input_output/ensight_output/EnSightOutput/GetEntityType.md)
- [VtkOutput](../../../../functions/input_output/vtk_output/VtkOutput/VtkOutput.md)
- [GetEntityType](../../../../functions/input_output/vtk_output/VtkOutput/GetEntityType.md)
- [PrintOutput](../../../../functions/input_output/vtk_output/VtkOutput/PrintOutput.md)
- [KRATOS_TEST_CASE_IN_SUITE](../../../../functions/mpi/tests/cpp_tests/sources/test_mpi_communicator/KRATOS_TEST_CASE_IN_SUITE.md)
- [__execute_test](../../../../functions/mpi/tests/test_distributed_import_model_part_utility/TestDistributedImportModelPartUtility/__execute_test.md)
- [test_with_real_model](../../../../functions/mpi/tests/test_distributed_model_part_initializer/TestDistributedModelPartInitializer/test_with_real_model.md)
- [test_GlobalNumberOf_Methods](../../../../functions/mpi/tests/test_mpi_communicator/TestMPICommunicator/test_GlobalNumberOf_Methods.md)
- [test_CombineModelParts](../../../../functions/tests/test_combine_model_part_modeler/TestCombineModelPartModeler/test_CombineModelParts.md)