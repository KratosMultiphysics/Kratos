---
type: C++ Method
title: ImportApplication
resource: sources/kernel.cpp#L151-L158
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/sources/kernel/Kernel/IsImported
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/python_interface/_ImportApplicationAsModule
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_interface/_ImportApplication
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kernel/Kernel/Initialize
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_utilities/test_suite/ImportApplicationIntoKernel
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void Kernel::ImportApplication(KratosApplication::Pointer pNewApplication)`

# Calls

- [IsImported](../../../../functions/sources/kernel/Kernel/IsImported.md)

# Called by

- [_ImportApplicationAsModule](../../../../functions/python_interface/_ImportApplicationAsModule.md)
- [_ImportApplication](../../../../functions/python_interface/_ImportApplication.md)
- [Initialize](../../../../functions/sources/kernel/Kernel/Initialize.md)
- [ImportApplicationIntoKernel](../../../../functions/tests/test_utilities/test_suite/ImportApplicationIntoKernel.md)