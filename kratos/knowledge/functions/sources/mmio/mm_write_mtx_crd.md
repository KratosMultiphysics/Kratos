---
type: C++ Function
title: mm_write_mtx_crd
resource: sources/mmio.c#L400-L441
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/sources/mmio/mm_typecode_to_str
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`int mm_write_mtx_crd(char fname[], int M, int N, int nz, int I[], int J[], double val[], MM_typecode matcode)`

# Calls

- [mm_typecode_to_str](../../../functions/sources/mmio/mm_typecode_to_str.md)