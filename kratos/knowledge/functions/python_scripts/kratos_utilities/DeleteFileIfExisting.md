---
type: Python Function
title: DeleteFileIfExisting
resource: python_scripts/kratos_utilities.py#L8-L15
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/sources/kratos_filesystem/remove
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/mpi/future/tests/test_mpi_vtu_output/check_pvtu
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/mpi/future/tests/test_mpi_vtu_output/check
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/mpi/future/tests/test_mpi_vtu_output/check_gauss
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/mpi/tests/test_mpi_serializer/TestMPISerializer/tearDown
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/compare_two_files_check_process/CompareTwoFilesCheckProcess/ExecuteFinalize
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/kratos_utilities/DeleteFilesEndingWith
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/restart_utility/RestartUtility/_ClearObsoleteRestartFiles
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_calculate_nodal_distance_to_skin_process/RemoveFiles
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_combine_model_part_modeler/TestCombineModelPartModeler/tearDown
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_compute_nodal_gradient_process/RemoveFiles
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_constraint_restart/TestConstraintRestart/tearDown
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_coordinate_transformation_utils/RemoveFiles
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_file_logger/TestFileLogger/tearDownClass
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_geometrical_object_bins/RemoveFiles
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_gid_io/TestGidIO/tearDown
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_gid_io_gauss_points/TestGiDIOGaussPoints/tearDown
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_import_model_part_modeler/TestImportMDPAModeler/tearDown
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_json_output_process/TestJsonOutputProcess/tearDown
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_kratos_parameters/TestParameters/test_file_serialization
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_matrix_market_interface/TestMatrixMarketInterface/test_mm_matrix_io
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_matrix_market_interface/TestMatrixMarketInterface/test_mm_vector_io
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_matrix_market_interface/TestMatrixMarketInterface/test_mm_matrix_io_cplx
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_matrix_market_interface/TestMatrixMarketInterface/test_mm_vector_io_cplx
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_model/TestModel/test_model_serialization
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_model_part_io/TestModelPartIO/tearDown
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_model_part_operation_utilities/TestModelPartOperationUtilities/tearDownClass
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_normal_utils/RemoveFiles
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_obj_io/RemoveFiles
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_restart/TestRestart/tearDown
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_sensitivity_utilities/RemoveFiles
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_serializer/TestSerializer/test_DataOnly
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_skin_detection_process/RemoveFiles
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_spatial_search/RemoveFiles
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_stl_io/RemoveFiles
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_table_logger/TestTableLogger/tearDown
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_vtu_output/TestVtuOutputBase/test_WriteMeshAsciiWithError
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_vtu_output/TestVtuOutputBase/test_WriteMeshBinaryWithError
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_vtu_output/TestVtuOutputBase/test_WriteMeshRawWithError
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_vtu_output/TestVtuOutputBase/test_WriteMeshCompressedRawWithError
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_vtu_output/TestVtuOutput/CheckVtuFile
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_vtu_output/TestVtuOutput/CheckVtmFile
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_vtu_output/TestVtuOutput/CheckSeriesFile
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_vtu_output/TestVtuOutput/CheckGaussVtuFile
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`def DeleteFileIfExisting(file_name):`

# Calls

- [remove](../../../functions/sources/kratos_filesystem/remove.md)

# Called by

- [check_pvtu](../../../functions/mpi/future/tests/test_mpi_vtu_output/check_pvtu.md)
- [check](../../../functions/mpi/future/tests/test_mpi_vtu_output/check.md)
- [check_gauss](../../../functions/mpi/future/tests/test_mpi_vtu_output/check_gauss.md)
- [tearDown](../../../functions/mpi/tests/test_mpi_serializer/TestMPISerializer/tearDown.md)
- [ExecuteFinalize](../../../functions/python_scripts/compare_two_files_check_process/CompareTwoFilesCheckProcess/ExecuteFinalize.md)
- [DeleteFilesEndingWith](../../../functions/python_scripts/kratos_utilities/DeleteFilesEndingWith.md)
- [_ClearObsoleteRestartFiles](../../../functions/python_scripts/restart_utility/RestartUtility/_ClearObsoleteRestartFiles.md)
- [RemoveFiles](../../../functions/tests/test_calculate_nodal_distance_to_skin_process/RemoveFiles.md)
- [tearDown](../../../functions/tests/test_combine_model_part_modeler/TestCombineModelPartModeler/tearDown.md)
- [RemoveFiles](../../../functions/tests/test_compute_nodal_gradient_process/RemoveFiles.md)
- [tearDown](../../../functions/tests/test_constraint_restart/TestConstraintRestart/tearDown.md)
- [RemoveFiles](../../../functions/tests/test_coordinate_transformation_utils/RemoveFiles.md)
- [tearDownClass](../../../functions/tests/test_file_logger/TestFileLogger/tearDownClass.md)
- [RemoveFiles](../../../functions/tests/test_geometrical_object_bins/RemoveFiles.md)
- [tearDown](../../../functions/tests/test_gid_io/TestGidIO/tearDown.md)
- [tearDown](../../../functions/tests/test_gid_io_gauss_points/TestGiDIOGaussPoints/tearDown.md)
- [tearDown](../../../functions/tests/test_import_model_part_modeler/TestImportMDPAModeler/tearDown.md)
- [tearDown](../../../functions/tests/test_json_output_process/TestJsonOutputProcess/tearDown.md)
- [test_file_serialization](../../../functions/tests/test_kratos_parameters/TestParameters/test_file_serialization.md)
- [test_mm_matrix_io](../../../functions/tests/test_matrix_market_interface/TestMatrixMarketInterface/test_mm_matrix_io.md)
- [test_mm_vector_io](../../../functions/tests/test_matrix_market_interface/TestMatrixMarketInterface/test_mm_vector_io.md)
- [test_mm_matrix_io_cplx](../../../functions/tests/test_matrix_market_interface/TestMatrixMarketInterface/test_mm_matrix_io_cplx.md)
- [test_mm_vector_io_cplx](../../../functions/tests/test_matrix_market_interface/TestMatrixMarketInterface/test_mm_vector_io_cplx.md)
- [test_model_serialization](../../../functions/tests/test_model/TestModel/test_model_serialization.md)
- [tearDown](../../../functions/tests/test_model_part_io/TestModelPartIO/tearDown.md)
- [tearDownClass](../../../functions/tests/test_model_part_operation_utilities/TestModelPartOperationUtilities/tearDownClass.md)
- [RemoveFiles](../../../functions/tests/test_normal_utils/RemoveFiles.md)
- [RemoveFiles](../../../functions/tests/test_obj_io/RemoveFiles.md)
- [tearDown](../../../functions/tests/test_restart/TestRestart/tearDown.md)
- [RemoveFiles](../../../functions/tests/test_sensitivity_utilities/RemoveFiles.md)
- [test_DataOnly](../../../functions/tests/test_serializer/TestSerializer/test_DataOnly.md)
- [RemoveFiles](../../../functions/tests/test_skin_detection_process/RemoveFiles.md)
- [RemoveFiles](../../../functions/tests/test_spatial_search/RemoveFiles.md)
- [RemoveFiles](../../../functions/tests/test_stl_io/RemoveFiles.md)
- [tearDown](../../../functions/tests/test_table_logger/TestTableLogger/tearDown.md)
- [test_WriteMeshAsciiWithError](../../../functions/tests/test_vtu_output/TestVtuOutputBase/test_WriteMeshAsciiWithError.md)
- [test_WriteMeshBinaryWithError](../../../functions/tests/test_vtu_output/TestVtuOutputBase/test_WriteMeshBinaryWithError.md)
- [test_WriteMeshRawWithError](../../../functions/tests/test_vtu_output/TestVtuOutputBase/test_WriteMeshRawWithError.md)
- [test_WriteMeshCompressedRawWithError](../../../functions/tests/test_vtu_output/TestVtuOutputBase/test_WriteMeshCompressedRawWithError.md)
- [CheckVtuFile](../../../functions/tests/test_vtu_output/TestVtuOutput/CheckVtuFile.md)
- [CheckVtmFile](../../../functions/tests/test_vtu_output/TestVtuOutput/CheckVtmFile.md)
- [CheckSeriesFile](../../../functions/tests/test_vtu_output/TestVtuOutput/CheckSeriesFile.md)
- [CheckGaussVtuFile](../../../functions/tests/test_vtu_output/TestVtuOutput/CheckGaussVtuFile.md)