---
type: Python Method
title: assertVectorAlmostEqual
resource: python_scripts/KratosUnittest.py#L77-L93
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/mpi/tests/test_distributed_sparse_matrices/TestDistributedSparseMatrices/test_matrix_construction
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/mpi/tests/test_mpi_data_communicator_python/TestMPIDataCommunicatorPython/testBroadcastOperations
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/mpi/tests/test_mpi_data_communicator_python/check_value_in_rank
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/mpi/tests/test_mpi_processes/TestMPIProcesses/testDiscontinuousDistanceProcessCutOnEdge2D
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/mpi/tests/test_mpi_processes/TestMPIProcesses/testDiscontinuousDistanceProcessCutThroughNode2D
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_json_output_process/TestJsonOutputProcess/test_output_of_vector_variable
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_json_output_process/TestJsonOutputProcess/test_output_of_resultant_vector_variable
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_kratos_parameters/TestParameters/test_vector_interface
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_matrix_interface/TestMatrixInterface/test_copy
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_matrix_market_interface/TestMatrixMarketInterface/test_mm_vector_io
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_matrix_market_interface/TestMatrixMarketInterface/test_mm_vector_io_cplx
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_model_part_io/TestModelPartIOMPI/test_model_part_io_read_entity_data
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_model_part_io/CheckEntities
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_normal_utils/check_normal_sensitivities
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_processes/TestProcesses/test_assign_master_slave_constraints_to_neighbours_process
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_rbf_shape_functions_utility/TestRBFShapeFunctionsUtility/test_rbf_shape_functions_utility_ill_conditioned
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_serializer/comparator
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_sparse_matrices/TestSparseMatrixInterface/test_assemble_with_csr_indices
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_sparse_matrices/TestSparseMatrixInterface/test_assemble_with_equation_ids
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_sparse_matrices/TestSparseMatrixInterface/test_assemble_with_csr_indices_int64
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_sparse_matrices/TestSparseMatrixInterface/test_assemble_with_equation_ids_int64
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_tensor_adaptors/TestTensorAdaptors/test_HistoricalVariableTensorAdaptorClone
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_tensor_adaptors/TestTensorAdaptors/test_VariableTensorAdaptorClone
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_tensor_adaptors/TestTensorAdaptors/test_NodePositionTensorAdaptorClone
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_tensor_adaptors/TestTensorAdaptors/test_CombinedTensorAdaptorClone
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_tensor_adaptors/TestTensorAdaptors/test_NodePositionTensorAdaptor1
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_tensor_adaptors/TestTensorAdaptors/test_NodePositionTensorAdaptor2
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_tensor_adaptors/TestTensorAdaptors/test_CombinedTensorAdaptorAxis1
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_tensor_adaptors/TestTensorAdaptors/test_CombinedTensorAdaptorAxis0
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_tensor_adaptors/TestTensorAdaptors/test_CombinedTensorAdaptorRavel
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_tensor_adaptors/TestTensorAdaptors/__TestEquationIdsTensorAdaptor
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_tensor_adaptors/TestTensorAdaptors/__TestTensorAdaptorShape
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_tensor_adaptors/TestTensorAdaptors/__TestCustomShape
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_tensor_adaptors/TestTensorAdaptors/__CheckValues
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_time_discretization/TestTimeDiscretization/_test_bdf_time_discretization
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_variable_utils/TestVariableUtils/test_copy_model_part_nodal_var
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_vector_interface/TestVectorInterface/test_sum
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`def assertVectorAlmostEqual(self, vector1, vector2, places=7, msg=None, delta=None):`

# Called by

- [test_matrix_construction](../../../../functions/mpi/tests/test_distributed_sparse_matrices/TestDistributedSparseMatrices/test_matrix_construction.md)
- [testBroadcastOperations](../../../../functions/mpi/tests/test_mpi_data_communicator_python/TestMPIDataCommunicatorPython/testBroadcastOperations.md)
- [check_value_in_rank](../../../../functions/mpi/tests/test_mpi_data_communicator_python/check_value_in_rank.md)
- [testDiscontinuousDistanceProcessCutOnEdge2D](../../../../functions/mpi/tests/test_mpi_processes/TestMPIProcesses/testDiscontinuousDistanceProcessCutOnEdge2D.md)
- [testDiscontinuousDistanceProcessCutThroughNode2D](../../../../functions/mpi/tests/test_mpi_processes/TestMPIProcesses/testDiscontinuousDistanceProcessCutThroughNode2D.md)
- [test_output_of_vector_variable](../../../../functions/tests/test_json_output_process/TestJsonOutputProcess/test_output_of_vector_variable.md)
- [test_output_of_resultant_vector_variable](../../../../functions/tests/test_json_output_process/TestJsonOutputProcess/test_output_of_resultant_vector_variable.md)
- [test_vector_interface](../../../../functions/tests/test_kratos_parameters/TestParameters/test_vector_interface.md)
- [test_copy](../../../../functions/tests/test_matrix_interface/TestMatrixInterface/test_copy.md)
- [test_mm_vector_io](../../../../functions/tests/test_matrix_market_interface/TestMatrixMarketInterface/test_mm_vector_io.md)
- [test_mm_vector_io_cplx](../../../../functions/tests/test_matrix_market_interface/TestMatrixMarketInterface/test_mm_vector_io_cplx.md)
- [test_model_part_io_read_entity_data](../../../../functions/tests/test_model_part_io/TestModelPartIOMPI/test_model_part_io_read_entity_data.md)
- [CheckEntities](../../../../functions/tests/test_model_part_io/CheckEntities.md)
- [check_normal_sensitivities](../../../../functions/tests/test_normal_utils/check_normal_sensitivities.md)
- [test_assign_master_slave_constraints_to_neighbours_process](../../../../functions/tests/test_processes/TestProcesses/test_assign_master_slave_constraints_to_neighbours_process.md)
- [test_rbf_shape_functions_utility_ill_conditioned](../../../../functions/tests/test_rbf_shape_functions_utility/TestRBFShapeFunctionsUtility/test_rbf_shape_functions_utility_ill_conditioned.md)
- [comparator](../../../../functions/tests/test_serializer/comparator.md)
- [test_assemble_with_csr_indices](../../../../functions/tests/test_sparse_matrices/TestSparseMatrixInterface/test_assemble_with_csr_indices.md)
- [test_assemble_with_equation_ids](../../../../functions/tests/test_sparse_matrices/TestSparseMatrixInterface/test_assemble_with_equation_ids.md)
- [test_assemble_with_csr_indices_int64](../../../../functions/tests/test_sparse_matrices/TestSparseMatrixInterface/test_assemble_with_csr_indices_int64.md)
- [test_assemble_with_equation_ids_int64](../../../../functions/tests/test_sparse_matrices/TestSparseMatrixInterface/test_assemble_with_equation_ids_int64.md)
- [test_HistoricalVariableTensorAdaptorClone](../../../../functions/tests/test_tensor_adaptors/TestTensorAdaptors/test_HistoricalVariableTensorAdaptorClone.md)
- [test_VariableTensorAdaptorClone](../../../../functions/tests/test_tensor_adaptors/TestTensorAdaptors/test_VariableTensorAdaptorClone.md)
- [test_NodePositionTensorAdaptorClone](../../../../functions/tests/test_tensor_adaptors/TestTensorAdaptors/test_NodePositionTensorAdaptorClone.md)
- [test_CombinedTensorAdaptorClone](../../../../functions/tests/test_tensor_adaptors/TestTensorAdaptors/test_CombinedTensorAdaptorClone.md)
- [test_NodePositionTensorAdaptor1](../../../../functions/tests/test_tensor_adaptors/TestTensorAdaptors/test_NodePositionTensorAdaptor1.md)
- [test_NodePositionTensorAdaptor2](../../../../functions/tests/test_tensor_adaptors/TestTensorAdaptors/test_NodePositionTensorAdaptor2.md)
- [test_CombinedTensorAdaptorAxis1](../../../../functions/tests/test_tensor_adaptors/TestTensorAdaptors/test_CombinedTensorAdaptorAxis1.md)
- [test_CombinedTensorAdaptorAxis0](../../../../functions/tests/test_tensor_adaptors/TestTensorAdaptors/test_CombinedTensorAdaptorAxis0.md)
- [test_CombinedTensorAdaptorRavel](../../../../functions/tests/test_tensor_adaptors/TestTensorAdaptors/test_CombinedTensorAdaptorRavel.md)
- [__TestEquationIdsTensorAdaptor](../../../../functions/tests/test_tensor_adaptors/TestTensorAdaptors/__TestEquationIdsTensorAdaptor.md)
- [__TestTensorAdaptorShape](../../../../functions/tests/test_tensor_adaptors/TestTensorAdaptors/__TestTensorAdaptorShape.md)
- [__TestCustomShape](../../../../functions/tests/test_tensor_adaptors/TestTensorAdaptors/__TestCustomShape.md)
- [__CheckValues](../../../../functions/tests/test_tensor_adaptors/TestTensorAdaptors/__CheckValues.md)
- [_test_bdf_time_discretization](../../../../functions/tests/test_time_discretization/TestTimeDiscretization/_test_bdf_time_discretization.md)
- [test_copy_model_part_nodal_var](../../../../functions/tests/test_variable_utils/TestVariableUtils/test_copy_model_part_nodal_var.md)
- [test_sum](../../../../functions/tests/test_vector_interface/TestVectorInterface/test_sum.md)