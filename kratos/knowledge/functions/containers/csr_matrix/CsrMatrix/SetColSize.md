---
type: C++ Method
title: SetColSize
resource: containers/csr_matrix.h#L310-L313
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/containers/distributed_csr_matrix/DistributedCsrMatrix/DistributedCsrMatrix
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/mpi/utilities/amgcl_distributed_csr_conversion_utilities/ConvertToCsrMatrix
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/amgcl_csr_conversion_utilities/ConvertToCsrMatrix
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void SetColSize(IndexType Ncols)`

# Called by

- [DistributedCsrMatrix](../../../../functions/containers/distributed_csr_matrix/DistributedCsrMatrix/DistributedCsrMatrix.md)
- [ConvertToCsrMatrix](../../../../functions/mpi/utilities/amgcl_distributed_csr_conversion_utilities/ConvertToCsrMatrix.md)
- [ConvertToCsrMatrix](../../../../functions/utilities/amgcl_csr_conversion_utilities/ConvertToCsrMatrix.md)