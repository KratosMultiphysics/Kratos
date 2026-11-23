---
type: C++ Method
title: CombinedTensorAdaptor
resource: tensor_adaptors/combined_tensor_adaptor.cpp#L41-L115
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/tensor_adaptors/combined_tensor_adaptor/GetModifiedTensorShape
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`CombinedTensorAdaptor<TDataType>::CombinedTensorAdaptor( const TensorAdaptorVectorType& rTensorAdaptorVector, const unsigned int Axis, const bool PerformCollectDataRecursively, const bool PerformStoreDataRecursively, const bool Copy) : mPerformCollectDataRecursively(PerformCollectDataRecursively), mPerformStoreDataRecursively(PerformStoreDataRecursively), mAxis(Axis)`

# Calls

- [GetModifiedTensorShape](../../../../functions/tensor_adaptors/combined_tensor_adaptor/GetModifiedTensorShape.md)