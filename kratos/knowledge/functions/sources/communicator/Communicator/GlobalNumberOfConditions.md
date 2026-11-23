---
type: C++ Method
title: GlobalNumberOfConditions
resource: sources/communicator.cpp#L111-L114
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/containers/edge_based_data_structure/EdgeBasedDataStructure/CalculateEdgeDataValues
    resolved_by: tree-sitter
    confidence: exact
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
  - target: functions/tests/test_skin_detection_process/TestSkinDetectionProcess/test_SkinDetectionProcess
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_stl_io/TestStlIO/test_WriteStlIO
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`Communicator::SizeType Communicator::GlobalNumberOfConditions() const`

# Called by

- [CalculateEdgeDataValues](../../../../functions/containers/edge_based_data_structure/EdgeBasedDataStructure/CalculateEdgeDataValues.md)
- [EnSightOutput](../../../../functions/input_output/ensight_output/EnSightOutput/EnSightOutput.md)
- [GetEntityType](../../../../functions/input_output/ensight_output/EnSightOutput/GetEntityType.md)
- [VtkOutput](../../../../functions/input_output/vtk_output/VtkOutput/VtkOutput.md)
- [GetEntityType](../../../../functions/input_output/vtk_output/VtkOutput/GetEntityType.md)
- [PrintOutput](../../../../functions/input_output/vtk_output/VtkOutput/PrintOutput.md)
- [__execute_test](../../../../functions/mpi/tests/test_distributed_import_model_part_utility/TestDistributedImportModelPartUtility/__execute_test.md)
- [test_with_real_model](../../../../functions/mpi/tests/test_distributed_model_part_initializer/TestDistributedModelPartInitializer/test_with_real_model.md)
- [test_GlobalNumberOf_Methods](../../../../functions/mpi/tests/test_mpi_communicator/TestMPICommunicator/test_GlobalNumberOf_Methods.md)
- [test_CombineModelParts](../../../../functions/tests/test_combine_model_part_modeler/TestCombineModelPartModeler/test_CombineModelParts.md)
- [test_SkinDetectionProcess](../../../../functions/tests/test_skin_detection_process/TestSkinDetectionProcess/test_SkinDetectionProcess.md)
- [test_WriteStlIO](../../../../functions/tests/test_stl_io/TestStlIO/test_WriteStlIO.md)