---
type: C++ Function
title: ImportApplicationIntoKernel
resource: tests/test_utilities/test_suite.h#L49-L53
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/sources/kernel/Kernel/IsImported
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kernel/Kernel/ImportApplication
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void ImportApplicationIntoKernel(KratosApplication::Pointer pNewApplication)`

# Calls

- [IsImported](../../../../functions/sources/kernel/Kernel/IsImported.md)
- [ImportApplication](../../../../functions/sources/kernel/Kernel/ImportApplication.md)