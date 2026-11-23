---
type: C++ Function
title: Transpose
resource: mpi/utilities/amgcl_distributed_csr_conversion_utilities.h#L144-L157
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/containers/distributed_csr_matrix/DistributedCsrMatrix/GetOffDiagonalIndex2DataInGlobalNumbering
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`static typename DistributedCsrMatrix<TDataType, TIndexType>::Pointer Transpose( const DistributedCsrMatrix<TDataType, TIndexType>& rA )`

# Calls

- [GetOffDiagonalIndex2DataInGlobalNumbering](../../../../functions/containers/distributed_csr_matrix/DistributedCsrMatrix/GetOffDiagonalIndex2DataInGlobalNumbering.md)