---
type: C++ Function
title: mm_read_unsymmetric_sparse
resource: sources/mmio.c#L17-L85
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/sources/mmio/mm_read_banner
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/mmio/mm_typecode_to_str
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/mmio/mm_read_mtx_crd_size
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`int mm_read_unsymmetric_sparse(const char *fname, int *M_, int *N_, int *nz_, double **val_, int **I_, int **J_)`

# Calls

- [mm_read_banner](../../../functions/sources/mmio/mm_read_banner.md)
- [mm_typecode_to_str](../../../functions/sources/mmio/mm_typecode_to_str.md)
- [mm_read_mtx_crd_size](../../../functions/sources/mmio/mm_read_mtx_crd_size.md)