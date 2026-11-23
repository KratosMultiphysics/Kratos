---
type: C++ Method
title: StoreData
resource: tensor_adaptors/historical_variable_tensor_adaptor.cpp#L170-L213
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/python/variable_indexing_python/visit
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tensor_adaptors/tensor_adaptor_utils/GetZeroValue
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tensor_adaptors/tensor_adaptor/TensorAdaptor<TDataType>/DataShape
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void HistoricalVariableTensorAdaptor::StoreData()`

# Calls

- [visit](../../../../functions/python/variable_indexing_python/visit.md)
- [GetZeroValue](../../../../functions/tensor_adaptors/tensor_adaptor_utils/GetZeroValue.md)
- [DataShape](../../../../functions/tensor_adaptors/tensor_adaptor/TensorAdaptor<TDataType>/DataShape.md)