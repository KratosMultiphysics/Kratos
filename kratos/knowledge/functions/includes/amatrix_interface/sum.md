---
type: C++ Function
title: sum
resource: includes/amatrix_interface.h#L675-L684
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/mpi/tests/test_mpi_serializer/executeComputeArea_Task
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_geometries_tensor_adaptor/TestGeometriesTensorAdaptor/test_ShapeFunctions_Gauss2
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_geometries_tensor_adaptor/TestGeometriesTensorAdaptor/test_ShapeFunctionsDerivatives_Gauss1
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_geometries_tensor_adaptor/TestGeometriesTensorAdaptor/test_ShapeFunctionsDerivatives_Gauss2
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_mls_shape_functions_utility/TestMLSShapeFunctionsUtility/testMLSShapeFunctionsUtility1x1Square1stOrder
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_mls_shape_functions_utility/TestMLSShapeFunctionsUtility/testMLSShapeFunctionsUtility1x1Square2ndOrder
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_nd_data/TestNDData/test_NDData1
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_nd_data/TestNDData/test_NDData2
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_nd_data/TestNDData/test_NDData3
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_rbf_shape_functions_utility/TestRBFShapeFunctionsUtility/test_rbf_shape_functions_utility_1x1_square
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_tensor_adaptors/TestTensorAdaptors/test_TensorAdaptorWithNDData
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_vectorized_interpolation/TestVectorizedInterpolation/Interpolate
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_voxel_modeler/TestVoxelMeshModeler/test_TestVoxelMeshModelerInRadius
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/rbf_shape_functions_utility/RBFShapeFunctionsUtility/CalculateShapeFunctions
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`typename TExpressionType::data_type sum( AMatrix::MatrixExpression<TExpressionType, TCategory> const& TheExpression)`

# Called by

- [executeComputeArea_Task](../../../functions/mpi/tests/test_mpi_serializer/executeComputeArea_Task.md)
- [test_ShapeFunctions_Gauss2](../../../functions/tests/test_geometries_tensor_adaptor/TestGeometriesTensorAdaptor/test_ShapeFunctions_Gauss2.md)
- [test_ShapeFunctionsDerivatives_Gauss1](../../../functions/tests/test_geometries_tensor_adaptor/TestGeometriesTensorAdaptor/test_ShapeFunctionsDerivatives_Gauss1.md)
- [test_ShapeFunctionsDerivatives_Gauss2](../../../functions/tests/test_geometries_tensor_adaptor/TestGeometriesTensorAdaptor/test_ShapeFunctionsDerivatives_Gauss2.md)
- [testMLSShapeFunctionsUtility1x1Square1stOrder](../../../functions/tests/test_mls_shape_functions_utility/TestMLSShapeFunctionsUtility/testMLSShapeFunctionsUtility1x1Square1stOrder.md)
- [testMLSShapeFunctionsUtility1x1Square2ndOrder](../../../functions/tests/test_mls_shape_functions_utility/TestMLSShapeFunctionsUtility/testMLSShapeFunctionsUtility1x1Square2ndOrder.md)
- [test_NDData1](../../../functions/tests/test_nd_data/TestNDData/test_NDData1.md)
- [test_NDData2](../../../functions/tests/test_nd_data/TestNDData/test_NDData2.md)
- [test_NDData3](../../../functions/tests/test_nd_data/TestNDData/test_NDData3.md)
- [test_rbf_shape_functions_utility_1x1_square](../../../functions/tests/test_rbf_shape_functions_utility/TestRBFShapeFunctionsUtility/test_rbf_shape_functions_utility_1x1_square.md)
- [test_TensorAdaptorWithNDData](../../../functions/tests/test_tensor_adaptors/TestTensorAdaptors/test_TensorAdaptorWithNDData.md)
- [Interpolate](../../../functions/tests/test_vectorized_interpolation/TestVectorizedInterpolation/Interpolate.md)
- [test_TestVoxelMeshModelerInRadius](../../../functions/tests/test_voxel_modeler/TestVoxelMeshModeler/test_TestVoxelMeshModelerInRadius.md)
- [CalculateShapeFunctions](../../../functions/utilities/rbf_shape_functions_utility/RBFShapeFunctionsUtility/CalculateShapeFunctions.md)