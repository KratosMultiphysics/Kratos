---
type: C++ Method
title: DataShape
resource: tensor_adaptors/tensor_adaptor.cpp#L135-L141
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/tensor_adaptors/gauss_point_variable_tensor_adaptor/GaussPointVariableTensorAdaptor/GaussPointVariableTensorAdaptor
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tensor_adaptors/historical_variable_tensor_adaptor/HistoricalVariableTensorAdaptor/HistoricalVariableTensorAdaptor
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tensor_adaptors/historical_variable_tensor_adaptor/HistoricalVariableTensorAdaptor/StoreData
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tensor_adaptors/node_position_tensor_adaptor/NodePositionTensorAdaptor/NodePositionTensorAdaptor
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tensor_adaptors/variable_tensor_adaptor/VariableTensorAdaptor/VariableTensorAdaptor
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tensor_adaptors/variable_tensor_adaptor/VariableTensorAdaptor/StoreData
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_tensor_adaptors/TestTensorAdaptors/__TestVariableTensorAdaptor
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_tensor_adaptors/TestTensorAdaptors/__TestVariableTensorAdaptorPartial
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_tensor_adaptors/TestTensorAdaptors/__TestTensorAdaptorShape
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`DenseVector<unsigned int> TensorAdaptor<TDataType>::DataShape() const`

# Called by

- [GaussPointVariableTensorAdaptor](../../../../functions/tensor_adaptors/gauss_point_variable_tensor_adaptor/GaussPointVariableTensorAdaptor/GaussPointVariableTensorAdaptor.md)
- [HistoricalVariableTensorAdaptor](../../../../functions/tensor_adaptors/historical_variable_tensor_adaptor/HistoricalVariableTensorAdaptor/HistoricalVariableTensorAdaptor.md)
- [StoreData](../../../../functions/tensor_adaptors/historical_variable_tensor_adaptor/HistoricalVariableTensorAdaptor/StoreData.md)
- [NodePositionTensorAdaptor](../../../../functions/tensor_adaptors/node_position_tensor_adaptor/NodePositionTensorAdaptor/NodePositionTensorAdaptor.md)
- [VariableTensorAdaptor](../../../../functions/tensor_adaptors/variable_tensor_adaptor/VariableTensorAdaptor/VariableTensorAdaptor.md)
- [StoreData](../../../../functions/tensor_adaptors/variable_tensor_adaptor/VariableTensorAdaptor/StoreData.md)
- [__TestVariableTensorAdaptor](../../../../functions/tests/test_tensor_adaptors/TestTensorAdaptors/__TestVariableTensorAdaptor.md)
- [__TestVariableTensorAdaptorPartial](../../../../functions/tests/test_tensor_adaptors/TestTensorAdaptors/__TestVariableTensorAdaptorPartial.md)
- [__TestTensorAdaptorShape](../../../../functions/tests/test_tensor_adaptors/TestTensorAdaptors/__TestTensorAdaptorShape.md)