---
type: C++ Method
title: local_size1
resource: containers/distributed_csr_matrix.h#L314-L317
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/containers/distributed_csr_matrix/DistributedCsrMatrix/ApplyHomogeneousDirichlet
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/containers/distributed_csr_matrix/DistributedCsrMatrix/ToMap
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/mpi/utilities/amgcl_distributed_csr_conversion_utilities/ConvertToAmgcl
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/petsc_conversion_tools/PetscCSRWrapper/__init__
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`TIndexType local_size1() const`

# Called by

- [ApplyHomogeneousDirichlet](../../../../functions/containers/distributed_csr_matrix/DistributedCsrMatrix/ApplyHomogeneousDirichlet.md)
- [ToMap](../../../../functions/containers/distributed_csr_matrix/DistributedCsrMatrix/ToMap.md)
- [ConvertToAmgcl](../../../../functions/mpi/utilities/amgcl_distributed_csr_conversion_utilities/ConvertToAmgcl.md)
- [__init__](../../../../functions/python_scripts/petsc_conversion_tools/PetscCSRWrapper/__init__.md)