---
type: C++ Module
title: registry
resource: sources/registry.cpp#L1-L185
generated:
  by: okf-rs/0.7.0
relationships:
  imports:
  - target: external/mutex
    resolved_by: tree-sitter
    confidence: exact
  - target: external/includes-registry-h
    resolved_by: tree-sitter
    confidence: exact
  member_of:
  - target: packages/KratosMultiphysics
    resolved_by: tree-sitter
    confidence: exact
---

# Contains

- [RemoveItem](../../functions/sources/registry/Registry/RemoveItem.md)
- [SetCurrentSource](../../functions/sources/registry/Registry/SetCurrentSource.md)
- [GetCurrentSource](../../functions/sources/registry/Registry/GetCurrentSource.md)
- [size](../../functions/sources/registry/Registry/size.md)
- [HasItem](../../functions/sources/registry/Registry/HasItem.md)
- [HasValue](../../functions/sources/registry/Registry/HasValue.md)
- [HasItems](../../functions/sources/registry/Registry/HasItems.md)
- [Info](../../functions/sources/registry/Registry/Info.md)
- [PrintInfo](../../functions/sources/registry/Registry/PrintInfo.md)
- [PrintData](../../functions/sources/registry/Registry/PrintData.md)
- [ToJson](../../functions/sources/registry/Registry/ToJson.md)
- [NotFoundError](../../functions/sources/registry/Registry/NotFoundError.md)

# Imports

- `<mutex>`
- `"includes/registry.h"`

# Member of

- [KratosMultiphysics](../../packages/KratosMultiphysics.md)