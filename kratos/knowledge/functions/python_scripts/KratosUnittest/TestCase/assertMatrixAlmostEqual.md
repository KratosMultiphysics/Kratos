---
type: Python Method
title: assertMatrixAlmostEqual
resource: python_scripts/KratosUnittest.py#L95-L128
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/mpi/tests/test_mpi_data_communicator_python/TestMPIDataCommunicatorPython/testBroadcastOperations
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/mpi/tests/test_mpi_data_communicator_python/check_value_in_rank
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_coordinate_transformation_utils/CheckSensitivities
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_kratos_parameters/TestParameters/test_matrix_interface
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_matrix_interface/TestMatrixInterface/test_copy
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_matrix_interface/TestMatrixInterface/test_matrix_sum
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_matrix_interface/TestMatrixInterface/test_matrix_diff
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_matrix_interface/TestMatrixInterface/test_scalar_prod
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_matrix_interface/TestMatrixInterface/test_empty_list_construction
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_matrix_interface/TestMatrixInterface/test_empty_list_of_list_construction
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_matrix_interface/TestMatrixInterface/test_buffer_protocol
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_matrix_market_interface/TestMatrixMarketInterface/test_mm_matrix_io
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_matrix_market_interface/TestMatrixMarketInterface/test_mm_matrix_io_cplx
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_model_part_io/TestModelPartIOMPI/test_model_part_io_read_entity_data
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_model_part_io/CheckEntities
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_sensitivity_utilities/TestSensitivityUtilitiesTwoDimSymmetricalSquare/__CheckAssignEntityDerivativesToNodesValues
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_serializer/comparator
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_tensor_adaptors/TestTensorAdaptors/__CheckValues
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`def assertMatrixAlmostEqual(self, matrix1, matrix2, places=7, msg=None, delta=None):`

# Called by

- [testBroadcastOperations](../../../../functions/mpi/tests/test_mpi_data_communicator_python/TestMPIDataCommunicatorPython/testBroadcastOperations.md)
- [check_value_in_rank](../../../../functions/mpi/tests/test_mpi_data_communicator_python/check_value_in_rank.md)
- [CheckSensitivities](../../../../functions/tests/test_coordinate_transformation_utils/CheckSensitivities.md)
- [test_matrix_interface](../../../../functions/tests/test_kratos_parameters/TestParameters/test_matrix_interface.md)
- [test_copy](../../../../functions/tests/test_matrix_interface/TestMatrixInterface/test_copy.md)
- [test_matrix_sum](../../../../functions/tests/test_matrix_interface/TestMatrixInterface/test_matrix_sum.md)
- [test_matrix_diff](../../../../functions/tests/test_matrix_interface/TestMatrixInterface/test_matrix_diff.md)
- [test_scalar_prod](../../../../functions/tests/test_matrix_interface/TestMatrixInterface/test_scalar_prod.md)
- [test_empty_list_construction](../../../../functions/tests/test_matrix_interface/TestMatrixInterface/test_empty_list_construction.md)
- [test_empty_list_of_list_construction](../../../../functions/tests/test_matrix_interface/TestMatrixInterface/test_empty_list_of_list_construction.md)
- [test_buffer_protocol](../../../../functions/tests/test_matrix_interface/TestMatrixInterface/test_buffer_protocol.md)
- [test_mm_matrix_io](../../../../functions/tests/test_matrix_market_interface/TestMatrixMarketInterface/test_mm_matrix_io.md)
- [test_mm_matrix_io_cplx](../../../../functions/tests/test_matrix_market_interface/TestMatrixMarketInterface/test_mm_matrix_io_cplx.md)
- [test_model_part_io_read_entity_data](../../../../functions/tests/test_model_part_io/TestModelPartIOMPI/test_model_part_io_read_entity_data.md)
- [CheckEntities](../../../../functions/tests/test_model_part_io/CheckEntities.md)
- [__CheckAssignEntityDerivativesToNodesValues](../../../../functions/tests/test_sensitivity_utilities/TestSensitivityUtilitiesTwoDimSymmetricalSquare/__CheckAssignEntityDerivativesToNodesValues.md)
- [comparator](../../../../functions/tests/test_serializer/comparator.md)
- [__CheckValues](../../../../functions/tests/test_tensor_adaptors/TestTensorAdaptors/__CheckValues.md)