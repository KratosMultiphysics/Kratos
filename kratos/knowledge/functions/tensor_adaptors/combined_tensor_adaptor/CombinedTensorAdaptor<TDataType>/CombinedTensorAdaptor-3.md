---
type: C++ Method
title: CombinedTensorAdaptor
resource: tensor_adaptors/combined_tensor_adaptor.cpp#L152-L173
visibility: private
generated:
  by: okf-rs/0.7.0
---

# Signature

`CombinedTensorAdaptor<TDataType>::CombinedTensorAdaptor( const CombinedTensorAdaptor& rOther, const bool PerformCollectDataRecursively, const bool PerformStoreDataRecursively, const bool Copy) : BaseType(rOther, Copy), mPerformCollectDataRecursively(PerformCollectDataRecursively), mPerformStoreDataRecursively(PerformStoreDataRecursively), mAxis(rOther.mAxis)`