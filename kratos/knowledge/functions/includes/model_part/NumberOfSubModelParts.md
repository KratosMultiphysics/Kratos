---
type: C++ Function
title: NumberOfSubModelParts
resource: includes/model_part.h#L1633-L1636
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/containers/model/Model/GetModelPartNames
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
  - target: functions/sources/model_part/ModelPart/GetSubModelPartNames
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/model_part/ModelPart/PrintData
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/model_part/ModelPart/save
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/input_output/test_stl_io/KRATOS_TEST_CASE_IN_SUITE
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/sources/test_model_part/KRATOS_TEST_CASE_IN_SUITE
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/sources/test_model_part_io/KRATOS_TEST_CASE_IN_SUITE
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_import_model_part_modeler/TestImportMDPAModeler/testImportMDPAModeler
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_model_part/TestModelPart/test_model_part_sub_model_parts
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_model_part/TestModelPart/test_clear_model_part
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_model_part_io/TestModelPartIO/test_model_part_io_read_using_submodelpart
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_model_part_io/TestModelPartIO/execute_test_model_part_io_read_model_part
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_model_part_io/TestModelPartIO/test_model_part_io_read_model_part_mesh_only
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_restart/TestRestart/_check_modelpart
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/auxiliar_model_part_utilities/AuxiliarModelPartUtilities/CopySubModelPartStructure
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/auxiliar_model_part_utilities/AuxiliarModelPartUtilities/DeepCopySubModelPart
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/model_part_combination_utilities/ModelPartCombinationUtilities/RecursiveAddOfModelPartsToList
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/model_part_combination_utilities/ModelPartCombinationUtilities/RecursiveAddEntities
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`SizeType NumberOfSubModelParts() const`

# Called by

- [GetModelPartNames](../../../functions/containers/model/Model/GetModelPartNames.md)
- [KRATOS_TEST_CASE_IN_SUITE](../../../functions/mpi/tests/cpp_tests/utilities/test_distributed_model_part_initializer/KRATOS_TEST_CASE_IN_SUITE.md)
- [__execute_test](../../../functions/mpi/tests/test_distributed_import_model_part_utility/TestDistributedImportModelPartUtility/__execute_test.md)
- [test_with_real_model](../../../functions/mpi/tests/test_distributed_model_part_initializer/TestDistributedModelPartInitializer/test_with_real_model.md)
- [GetSubModelPartNames](../../../functions/sources/model_part/ModelPart/GetSubModelPartNames.md)
- [PrintData](../../../functions/sources/model_part/ModelPart/PrintData.md)
- [save](../../../functions/sources/model_part/ModelPart/save.md)
- [KRATOS_TEST_CASE_IN_SUITE](../../../functions/tests/cpp_tests/input_output/test_stl_io/KRATOS_TEST_CASE_IN_SUITE.md)
- [KRATOS_TEST_CASE_IN_SUITE](../../../functions/tests/cpp_tests/sources/test_model_part/KRATOS_TEST_CASE_IN_SUITE.md)
- [KRATOS_TEST_CASE_IN_SUITE](../../../functions/tests/cpp_tests/sources/test_model_part_io/KRATOS_TEST_CASE_IN_SUITE.md)
- [testImportMDPAModeler](../../../functions/tests/test_import_model_part_modeler/TestImportMDPAModeler/testImportMDPAModeler.md)
- [test_model_part_sub_model_parts](../../../functions/tests/test_model_part/TestModelPart/test_model_part_sub_model_parts.md)
- [test_clear_model_part](../../../functions/tests/test_model_part/TestModelPart/test_clear_model_part.md)
- [test_model_part_io_read_using_submodelpart](../../../functions/tests/test_model_part_io/TestModelPartIO/test_model_part_io_read_using_submodelpart.md)
- [execute_test_model_part_io_read_model_part](../../../functions/tests/test_model_part_io/TestModelPartIO/execute_test_model_part_io_read_model_part.md)
- [test_model_part_io_read_model_part_mesh_only](../../../functions/tests/test_model_part_io/TestModelPartIO/test_model_part_io_read_model_part_mesh_only.md)
- [_check_modelpart](../../../functions/tests/test_restart/TestRestart/_check_modelpart.md)
- [CopySubModelPartStructure](../../../functions/utilities/auxiliar_model_part_utilities/AuxiliarModelPartUtilities/CopySubModelPartStructure.md)
- [DeepCopySubModelPart](../../../functions/utilities/auxiliar_model_part_utilities/AuxiliarModelPartUtilities/DeepCopySubModelPart.md)
- [RecursiveAddOfModelPartsToList](../../../functions/utilities/model_part_combination_utilities/ModelPartCombinationUtilities/RecursiveAddOfModelPartsToList.md)
- [RecursiveAddEntities](../../../functions/utilities/model_part_combination_utilities/ModelPartCombinationUtilities/RecursiveAddEntities.md)