---
type: C++ Method
title: GetLibraryList
resource: sources/kernel.cpp#L118-L127
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/sources/kernel/Kernel/IsLibraryAvailable
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`std::unordered_set<std::string> Kernel::GetLibraryList()`

# Called by

- [IsLibraryAvailable](../../../../functions/sources/kernel/Kernel/IsLibraryAvailable.md)