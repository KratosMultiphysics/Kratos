---
type: C++ Function
title: SparseMultiply
resource: mpi/utilities/amgcl_distributed_csr_spmm_utilities.h#L41-L62
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/containers/distributed_csr_matrix/DistributedCsrMatrix/GetOffDiagonalIndex2DataInGlobalNumbering
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`static typename DistributedCsrMatrix<TDataType, TIndexType>::Pointer SparseMultiply( const DistributedCsrMatrix<TDataType, TIndexType>& rA, const DistributedCsrMatrix<TDataType, TIndexType>& rB )`

# Calls

- [GetOffDiagonalIndex2DataInGlobalNumbering](../../../../functions/containers/distributed_csr_matrix/DistributedCsrMatrix/GetOffDiagonalIndex2DataInGlobalNumbering.md)