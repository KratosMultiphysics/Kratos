---
type: C++ Method
title: IsLibraryAvailable
resource: sources/kernel.cpp#L133-L145
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/sources/kernel/Kernel/GetLibraryList
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/tests/test_kernel/TestKernel/testIsLibraryAvailable
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`bool Kernel::IsLibraryAvailable(const std::string& rLibraryName) const`

# Calls

- [GetLibraryList](../../../../functions/sources/kernel/Kernel/GetLibraryList.md)

# Called by

- [testIsLibraryAvailable](../../../../functions/tests/test_kernel/TestKernel/testIsLibraryAvailable.md)