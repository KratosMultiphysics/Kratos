---
type: C++ Function
title: WriteMatrixMarketMatrix
resource: sources/matrix_market_interface.cpp#L325-L437
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/sources/mmio/mm_write_banner
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/containers/csr_matrix/CsrMatrix/nnz
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/mmio/mm_write_mtx_crd_size
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`bool WriteMatrixMarketMatrix(const char *FileName, const CompressedMatrixType &M, bool Symmetric)`

# Calls

- [mm_write_banner](../../../functions/sources/mmio/mm_write_banner.md)
- [nnz](../../../functions/containers/csr_matrix/CsrMatrix/nnz.md)
- [mm_write_mtx_crd_size](../../../functions/sources/mmio/mm_write_mtx_crd_size.md)