---
type: C++ Method
title: TransposeSpMV
resource: containers/distributed_csr_matrix.h#L515-L533
generated:
  by: okf-rs/0.7.0
---

# Signature

`DistributedVectorExporter<TIndexType>* TransposeSpMV(const DistributedSystemVector<TDataType,TIndexType>& x, DistributedSystemVector<TDataType,TIndexType>& y, DistributedVectorExporter<TIndexType>* pTransposeExporter = nullptr ) const`