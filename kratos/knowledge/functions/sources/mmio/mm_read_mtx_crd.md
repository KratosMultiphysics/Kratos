---
type: C++ Function
title: mm_read_mtx_crd
resource: sources/mmio.c#L334-L385
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/sources/mmio/mm_read_banner
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/mmio/mm_is_valid
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/mmio/mm_read_mtx_crd_size
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/mmio/mm_read_mtx_crd_data
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`int mm_read_mtx_crd(char *fname, int *M, int *N, int *nz, int **I, int **J, double **val, MM_typecode *matcode)`

# Calls

- [mm_read_banner](../../../functions/sources/mmio/mm_read_banner.md)
- [mm_is_valid](../../../functions/sources/mmio/mm_is_valid.md)
- [mm_read_mtx_crd_size](../../../functions/sources/mmio/mm_read_mtx_crd_size.md)
- [mm_read_mtx_crd_data](../../../functions/sources/mmio/mm_read_mtx_crd_data.md)