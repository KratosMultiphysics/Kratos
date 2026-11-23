---
type: C++ Function
title: GetZeroValue
resource: tensor_adaptors/tensor_adaptor_utils.h#L239-L258
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/tensor_adaptors/historical_variable_tensor_adaptor/HistoricalVariableTensorAdaptor/StoreData
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tensor_adaptors/variable_tensor_adaptor/VariableTensorAdaptor/StoreData
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`[[nodiscard]] static TDataType GetZeroValue( const Variable<TDataType>& rVariable, const DenseVector<unsigned int>& rDataShape)`

# Called by

- [StoreData](../../../functions/tensor_adaptors/historical_variable_tensor_adaptor/HistoricalVariableTensorAdaptor/StoreData.md)
- [StoreData](../../../functions/tensor_adaptors/variable_tensor_adaptor/VariableTensorAdaptor/StoreData.md)