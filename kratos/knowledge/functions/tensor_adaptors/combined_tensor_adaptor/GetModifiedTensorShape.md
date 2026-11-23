---
type: C++ Function
title: GetModifiedTensorShape
resource: tensor_adaptors/combined_tensor_adaptor.cpp#L28-L37
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/tensor_adaptors/combined_tensor_adaptor/CombinedTensorAdaptor<TDataType>/CombinedTensorAdaptor
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tensor_adaptors/combined_tensor_adaptor/CombinedTensorAdaptor<TDataType>/CollectData
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`[[nodiscard]] DenseVector<unsigned int> GetModifiedTensorShape( const IndexType NumberOfDimensions, const DenseVector<unsigned int>& rCurrentTensorShape)`

# Called by

- [CombinedTensorAdaptor](../../../functions/tensor_adaptors/combined_tensor_adaptor/CombinedTensorAdaptor<TDataType>/CombinedTensorAdaptor.md)
- [CollectData](../../../functions/tensor_adaptors/combined_tensor_adaptor/CombinedTensorAdaptor<TDataType>/CollectData.md)