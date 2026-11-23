---
type: C++ Method
title: CombinedTensorAdaptor
resource: tensor_adaptors/combined_tensor_adaptor.cpp#L118-L149
visibility: private
generated:
  by: okf-rs/0.7.0
---

# Signature

`CombinedTensorAdaptor<TDataType>::CombinedTensorAdaptor( const TensorAdaptorVectorType& rTensorAdaptorVector, const bool PerformCollectDataRecursively, const bool PerformStoreDataRecursively, const bool Copy) : mPerformCollectDataRecursively(PerformCollectDataRecursively), mPerformStoreDataRecursively(PerformStoreDataRecursively), mAxis(-1) // represent that this combined tensor adaptor is used with ravel`