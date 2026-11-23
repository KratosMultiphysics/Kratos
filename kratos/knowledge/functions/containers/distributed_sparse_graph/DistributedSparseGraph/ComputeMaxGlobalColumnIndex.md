---
type: C++ Method
title: ComputeMaxGlobalColumnIndex
resource: containers/distributed_sparse_graph.h#L146-L151
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/containers/distributed_sparse_graph/DistributedSparseGraph/ComputeLocalMinMaxColumnIndex
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/containers/distributed_csr_matrix/DistributedCsrMatrix/DistributedCsrMatrix
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`IndexType ComputeMaxGlobalColumnIndex() const`

# Calls

- [ComputeLocalMinMaxColumnIndex](../../../../functions/containers/distributed_sparse_graph/DistributedSparseGraph/ComputeLocalMinMaxColumnIndex.md)

# Called by

- [DistributedCsrMatrix](../../../../functions/containers/distributed_csr_matrix/DistributedCsrMatrix/DistributedCsrMatrix.md)