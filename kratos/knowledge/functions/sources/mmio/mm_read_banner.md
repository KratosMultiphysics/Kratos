---
type: C++ Function
title: mm_read_banner
resource: sources/mmio.c#L97-L180
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/sources/matrix_market_interface/ReadMatrixMarketMatrix
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/matrix_market_interface/ReadMatrixMarketVector
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/mmio/mm_read_unsymmetric_sparse
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/mmio/mm_read_mtx_crd
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`int mm_read_banner(FILE *f, MM_typecode *matcode)`

# Called by

- [ReadMatrixMarketMatrix](../../../functions/sources/matrix_market_interface/ReadMatrixMarketMatrix.md)
- [ReadMatrixMarketVector](../../../functions/sources/matrix_market_interface/ReadMatrixMarketVector.md)
- [mm_read_unsymmetric_sparse](../../../functions/sources/mmio/mm_read_unsymmetric_sparse.md)
- [mm_read_mtx_crd](../../../functions/sources/mmio/mm_read_mtx_crd.md)