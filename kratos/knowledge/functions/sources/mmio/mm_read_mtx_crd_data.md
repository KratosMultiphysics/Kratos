---
type: C++ Function
title: mm_read_mtx_crd_data
resource: sources/mmio.c#L266-L297
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/sources/mmio/mm_read_mtx_crd
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`int mm_read_mtx_crd_data(FILE *f, int M, int N, int nz, int I[], int J[], double val[], MM_typecode matcode)`

# Called by

- [mm_read_mtx_crd](../../../functions/sources/mmio/mm_read_mtx_crd.md)