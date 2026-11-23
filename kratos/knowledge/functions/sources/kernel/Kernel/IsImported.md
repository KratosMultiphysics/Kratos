---
type: C++ Method
title: IsImported
resource: sources/kernel.cpp#L129-L131
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/sources/kernel/Kernel/Initialize
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kernel/Kernel/ImportApplication
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_importing/TestImporting/test_has_application
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_utilities/test_suite/ImportApplicationIntoKernel
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`bool Kernel::IsImported(const std::string& rApplicationName) const`

# Called by

- [Initialize](../../../../functions/sources/kernel/Kernel/Initialize.md)
- [ImportApplication](../../../../functions/sources/kernel/Kernel/ImportApplication.md)
- [test_has_application](../../../../functions/tests/test_importing/TestImporting/test_has_application.md)
- [ImportApplicationIntoKernel](../../../../functions/tests/test_utilities/test_suite/ImportApplicationIntoKernel.md)