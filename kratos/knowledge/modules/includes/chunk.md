---
type: C++ Module
title: chunk
resource: includes/chunk.h#L1-L338
generated:
  by: okf-rs/0.7.0
relationships:
  imports:
  - target: external/atomic
    resolved_by: tree-sitter
    confidence: exact
  - target: external/includes-checks-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/includes-lock-object-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/utilities-openmp-utils-h
    resolved_by: tree-sitter
    confidence: exact
  member_of:
  - target: packages/KratosMultiphysics
    resolved_by: tree-sitter
    confidence: exact
---

# Contains

- [Chunk](../../classes/includes/chunk/Chunk.md)
- [Chunk](../../functions/includes/chunk/Chunk/Chunk.md)
- [Chunk](../../functions/includes/chunk/Chunk/Chunk-2.md)
- [Chunk](../../functions/includes/chunk/Chunk/Chunk-3.md)
- [Chunk](../../functions/includes/chunk/Chunk/Chunk-4.md)
- [Initialize](../../functions/includes/chunk/Chunk/Initialize.md)
- [Allocate](../../functions/includes/chunk/Chunk/Allocate.md)
- [Deallocate](../../functions/includes/chunk/Chunk/Deallocate.md)
- [Release](../../functions/includes/chunk/Chunk/Release.md)
- [Size](../../functions/includes/chunk/Chunk/Size.md)
- [MemoryUsed](../../functions/includes/chunk/Chunk/MemoryUsed.md)
- [MemoryOverhead](../../functions/includes/chunk/Chunk/MemoryOverhead.md)
- [GetNumberOfBlocks](../../functions/includes/chunk/Chunk/GetNumberOfBlocks.md)
- [GetNumberOfAvailableBlocks](../../functions/includes/chunk/Chunk/GetNumberOfAvailableBlocks.md)
- [pGetData](../../functions/includes/chunk/Chunk/pGetData.md)
- [pDataBegin](../../functions/includes/chunk/Chunk/pDataBegin.md)
- [pDataEnd](../../functions/includes/chunk/Chunk/pDataEnd.md)
- [HasAvailableBlock](../../functions/includes/chunk/Chunk/HasAvailableBlock.md)
- [Has](../../functions/includes/chunk/Chunk/Has.md)
- [IsEmpty](../../functions/includes/chunk/Chunk/IsEmpty.md)
- [IsFull](../../functions/includes/chunk/Chunk/IsFull.md)
- [IsInitialized](../../functions/includes/chunk/Chunk/IsInitialized.md)
- [IsReleased](../../functions/includes/chunk/Chunk/IsReleased.md)
- [lock](../../functions/includes/chunk/Chunk/lock.md)
- [unlock](../../functions/includes/chunk/Chunk/unlock.md)
- [Info](../../functions/includes/chunk/Chunk/Info.md)
- [PrintInfo](../../functions/includes/chunk/Chunk/PrintInfo.md)
- [PrintData](../../functions/includes/chunk/Chunk/PrintData.md)
- [GetHeaderSize](../../functions/includes/chunk/Chunk/GetHeaderSize.md)
- [DataSize](../../functions/includes/chunk/Chunk/DataSize.md)
- [AllocatableDataSize](../../functions/includes/chunk/Chunk/AllocatableDataSize.md)
- [GetBlockSize](../../functions/includes/chunk/Chunk/GetBlockSize.md)

# Imports

- `<atomic>`
- `"includes/checks.h"`
- `"includes/lock_object.h"`
- `"utilities/openmp_utils.h"`

# Member of

- [KratosMultiphysics](../../packages/KratosMultiphysics.md)