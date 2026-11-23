---
type: C++ Method
title: HasSubModelPart
resource: sources/model_part.cpp#L2144-L2159
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/containers/model/Model/HasModelPart
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/modeler/combine_model_part_modeler/CombineModelPartModeler/DuplicateSubModelParts
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/modeler/connectivity_preserve_modeler/ConnectivityPreserveModeler/DuplicateSubModelParts
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/modeler/copy_properties_modeler/CopyPropertiesModeler/RecursivelyCopyProperties
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/modeler/duplicate_mesh_modeler/DuplicateMeshModeler/GenerateSubModelParts
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
  - target: functions/mpi/tests/test_mpi_communicator/TestMPICommunicator/_read_model_part_mpi
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/mpi/utilities/distributed_model_part_initializer/RecursiveCreateModelParts
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/split_internal_interfaces_process/SplitInternalInterfacesProcess
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/containers/test_model/KRATOS_TEST_CASE_IN_SUITE
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/input_output/test_stl_io/KRATOS_TEST_CASE_IN_SUITE
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/processes/test_structured_mesh_generator_process/KRATOS_TEST_CASE_IN_SUITE
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/sources/test_model_part/KRATOS_TEST_CASE_IN_SUITE
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/sources/test_model_part_io/KRATOS_TEST_CASE_IN_SUITE
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/utilities/test_auxiliar_model_part_utilities/KRATOS_TEST_CASE_IN_SUITE
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_duplicate_mesh_modeler/TestDuplicateMeshModeler/test_DuplicateMeshModeler
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_import_model_part_modeler/TestImportMDPAModeler/testImportMDPAModeler
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_model_part/TestModelPart/test_model_part_sub_model_parts
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_model_part_combination_utilities/TestModelPartCombinationUtilities/test_combine_model_parts
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
  - target: functions/utilities/shifted_boundary_meshless_interface_utility/ShiftedBoundaryMeshlessInterfaceUtility/ShiftedBoundaryMeshlessInterfaceUtility
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`bool ModelPart::HasSubModelPart(std::string const& ThisSubModelPartName) const`

# Called by

- [HasModelPart](../../../../functions/containers/model/Model/HasModelPart.md)
- [DuplicateSubModelParts](../../../../functions/modeler/combine_model_part_modeler/CombineModelPartModeler/DuplicateSubModelParts.md)
- [DuplicateSubModelParts](../../../../functions/modeler/connectivity_preserve_modeler/ConnectivityPreserveModeler/DuplicateSubModelParts.md)
- [RecursivelyCopyProperties](../../../../functions/modeler/copy_properties_modeler/CopyPropertiesModeler/RecursivelyCopyProperties.md)
- [GenerateSubModelParts](../../../../functions/modeler/duplicate_mesh_modeler/DuplicateMeshModeler/GenerateSubModelParts.md)
- [KRATOS_TEST_CASE_IN_SUITE](../../../../functions/mpi/tests/cpp_tests/utilities/test_distributed_model_part_initializer/KRATOS_TEST_CASE_IN_SUITE.md)
- [__execute_test](../../../../functions/mpi/tests/test_distributed_import_model_part_utility/TestDistributedImportModelPartUtility/__execute_test.md)
- [test_with_real_model](../../../../functions/mpi/tests/test_distributed_model_part_initializer/TestDistributedModelPartInitializer/test_with_real_model.md)
- [_read_model_part_mpi](../../../../functions/mpi/tests/test_mpi_communicator/TestMPICommunicator/_read_model_part_mpi.md)
- [RecursiveCreateModelParts](../../../../functions/mpi/utilities/distributed_model_part_initializer/RecursiveCreateModelParts.md)
- [SplitInternalInterfacesProcess](../../../../functions/processes/split_internal_interfaces_process/SplitInternalInterfacesProcess.md)
- [KRATOS_TEST_CASE_IN_SUITE](../../../../functions/tests/cpp_tests/containers/test_model/KRATOS_TEST_CASE_IN_SUITE.md)
- [KRATOS_TEST_CASE_IN_SUITE](../../../../functions/tests/cpp_tests/input_output/test_stl_io/KRATOS_TEST_CASE_IN_SUITE.md)
- [KRATOS_TEST_CASE_IN_SUITE](../../../../functions/tests/cpp_tests/processes/test_structured_mesh_generator_process/KRATOS_TEST_CASE_IN_SUITE.md)
- [KRATOS_TEST_CASE_IN_SUITE](../../../../functions/tests/cpp_tests/sources/test_model_part/KRATOS_TEST_CASE_IN_SUITE.md)
- [KRATOS_TEST_CASE_IN_SUITE](../../../../functions/tests/cpp_tests/sources/test_model_part_io/KRATOS_TEST_CASE_IN_SUITE.md)
- [KRATOS_TEST_CASE_IN_SUITE](../../../../functions/tests/cpp_tests/utilities/test_auxiliar_model_part_utilities/KRATOS_TEST_CASE_IN_SUITE.md)
- [test_DuplicateMeshModeler](../../../../functions/tests/test_duplicate_mesh_modeler/TestDuplicateMeshModeler/test_DuplicateMeshModeler.md)
- [testImportMDPAModeler](../../../../functions/tests/test_import_model_part_modeler/TestImportMDPAModeler/testImportMDPAModeler.md)
- [test_model_part_sub_model_parts](../../../../functions/tests/test_model_part/TestModelPart/test_model_part_sub_model_parts.md)
- [test_combine_model_parts](../../../../functions/tests/test_model_part_combination_utilities/TestModelPartCombinationUtilities/test_combine_model_parts.md)
- [execute_test_model_part_io_read_model_part](../../../../functions/tests/test_model_part_io/TestModelPartIO/execute_test_model_part_io_read_model_part.md)
- [test_model_part_io_read_model_part_mesh_only](../../../../functions/tests/test_model_part_io/TestModelPartIO/test_model_part_io_read_model_part_mesh_only.md)
- [_check_modelpart](../../../../functions/tests/test_restart/TestRestart/_check_modelpart.md)
- [ShiftedBoundaryMeshlessInterfaceUtility](../../../../functions/utilities/shifted_boundary_meshless_interface_utility/ShiftedBoundaryMeshlessInterfaceUtility/ShiftedBoundaryMeshlessInterfaceUtility.md)