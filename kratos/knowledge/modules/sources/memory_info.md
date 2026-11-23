---
type: C++ Module
title: memory_info
resource: sources/memory_info.cpp#L1-L116
generated:
  by: okf-rs/0.7.0
relationships:
  imports:
  - target: external/iostream
    resolved_by: tree-sitter
    confidence: exact
  - target: external/sstream
    resolved_by: tree-sitter
    confidence: exact
  - target: external/string
    resolved_by: tree-sitter
    confidence: exact
  - target: external/fstream
    resolved_by: tree-sitter
    confidence: exact
  - target: external/includes-memory-info-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/sys-resource-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/unistd-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/mach-mach-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/windows-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/psapi-h
    resolved_by: tree-sitter
    confidence: exact
  member_of:
  - target: packages/KratosMultiphysics
    resolved_by: tree-sitter
    confidence: exact
---

# Contains

- [GetPeakMemoryUsage](../../functions/sources/memory_info/MemoryInfo/GetPeakMemoryUsage.md)
- [GetCurrentMemoryUsage](../../functions/sources/memory_info/MemoryInfo/GetCurrentMemoryUsage.md)
- [HumanReadableSize](../../functions/sources/memory_info/MemoryInfo/HumanReadableSize.md)

# Imports

- `<iostream>`
- `<sstream>`
- `<string>`
- `<fstream>`
- `"includes/memory_info.h"`
- `<sys/resource.h>`
- `<unistd.h>`
- `<mach/mach.h>`
- `<windows.h>`
- `<psapi.h>`

# Member of

- [KratosMultiphysics](../../packages/KratosMultiphysics.md)