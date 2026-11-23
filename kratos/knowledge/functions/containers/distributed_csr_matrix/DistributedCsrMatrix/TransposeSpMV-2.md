---
type: C++ Method
title: TransposeSpMV
resource: containers/distributed_csr_matrix.h#L536-L557
generated:
  by: okf-rs/0.7.0
---

# Signature

`DistributedVectorExporter<TIndexType>* TransposeSpMV( TDataType alpha, const DistributedSystemVector<TDataType,TIndexType>& x, TDataType beta, DistributedSystemVector<TDataType,TIndexType>& y, DistributedVectorExporter<TIndexType>* pTransposeExporter = nullptr ) const`