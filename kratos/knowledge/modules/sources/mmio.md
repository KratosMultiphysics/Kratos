---
type: C++ Module
title: mmio
resource: sources/mmio.c#L1-L513
generated:
  by: okf-rs/0.7.0
relationships:
  imports:
  - target: external/stdio-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/string-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/stdlib-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/ctype-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/includes-mmio-h
    resolved_by: tree-sitter
    confidence: exact
  member_of:
  - target: packages/KratosMultiphysics
    resolved_by: tree-sitter
    confidence: exact
---

# Contains

- [mm_read_unsymmetric_sparse](../../functions/sources/mmio/mm_read_unsymmetric_sparse.md)
- [mm_is_valid](../../functions/sources/mmio/mm_is_valid.md)
- [mm_read_banner](../../functions/sources/mmio/mm_read_banner.md)
- [mm_write_mtx_crd_size](../../functions/sources/mmio/mm_write_mtx_crd_size.md)
- [mm_read_mtx_crd_size](../../functions/sources/mmio/mm_read_mtx_crd_size.md)
- [mm_read_mtx_array_size](../../functions/sources/mmio/mm_read_mtx_array_size.md)
- [mm_write_mtx_array_size](../../functions/sources/mmio/mm_write_mtx_array_size.md)
- [mm_read_mtx_crd_data](../../functions/sources/mmio/mm_read_mtx_crd_data.md)
- [mm_read_mtx_crd_entry](../../functions/sources/mmio/mm_read_mtx_crd_entry.md)
- [mm_read_mtx_crd](../../functions/sources/mmio/mm_read_mtx_crd.md)
- [mm_write_banner](../../functions/sources/mmio/mm_write_banner.md)
- [mm_write_mtx_crd](../../functions/sources/mmio/mm_write_mtx_crd.md)
- [mm_strdup](../../functions/sources/mmio/mm_strdup.md)
- [mm_typecode_to_str](../../functions/sources/mmio/mm_typecode_to_str.md)

# Imports

- `<stdio.h>`
- `<string.h>`
- `<stdlib.h>`
- `<ctype.h>`
- `"includes/mmio.h"`

# Member of

- [KratosMultiphysics](../../packages/KratosMultiphysics.md)