---
type: C++ Method
title: PrintData
resource: containers/distributed_csr_matrix.h#L823-L870
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/containers/csr_matrix/CsrMatrix/nnz
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/containers/distributed_csr_matrix/DistributedCsrMatrix/GetDiagonalIndex2DataInGlobalNumbering
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/containers/distributed_csr_matrix/DistributedCsrMatrix/GetOffDiagonalIndex2DataInGlobalNumbering
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void PrintData(std::ostream& rOStream) const`

# Calls

- [nnz](../../../../functions/containers/csr_matrix/CsrMatrix/nnz.md)
- [GetDiagonalIndex2DataInGlobalNumbering](../../../../functions/containers/distributed_csr_matrix/DistributedCsrMatrix/GetDiagonalIndex2DataInGlobalNumbering.md)
- [GetOffDiagonalIndex2DataInGlobalNumbering](../../../../functions/containers/distributed_csr_matrix/DistributedCsrMatrix/GetOffDiagonalIndex2DataInGlobalNumbering.md)