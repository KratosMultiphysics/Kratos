---
type: C++ Function
title: PrintData
resource: includes/memory_info.h#L91-L97
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/sources/memory_info/MemoryInfo/HumanReadableSize
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/memory_info/MemoryInfo/GetCurrentMemoryUsage
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/memory_info/MemoryInfo/GetPeakMemoryUsage
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`virtual void PrintData(std::ostream &rOStream) const`

# Calls

- [HumanReadableSize](../../../functions/sources/memory_info/MemoryInfo/HumanReadableSize.md)
- [GetCurrentMemoryUsage](../../../functions/sources/memory_info/MemoryInfo/GetCurrentMemoryUsage.md)
- [GetPeakMemoryUsage](../../../functions/sources/memory_info/MemoryInfo/GetPeakMemoryUsage.md)