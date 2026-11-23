---
type: Python Method
title: loadTestsFromTestCases
resource: python_scripts/KratosUnittest.py#L16-L26
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/future/tests/test_KratosFutureCore/AssembleTestSuites
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/mpi/future/tests/test_KratosFutureMPICore/AssembleTestSuites
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/mpi/tests/test_KratosMPICore/AssembleTestSuites
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_KratosCore/AssembleTestSuites
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`def loadTestsFromTestCases(self, testCaseClasses):`

# Called by

- [AssembleTestSuites](../../../../functions/future/tests/test_KratosFutureCore/AssembleTestSuites.md)
- [AssembleTestSuites](../../../../functions/mpi/future/tests/test_KratosFutureMPICore/AssembleTestSuites.md)
- [AssembleTestSuites](../../../../functions/mpi/tests/test_KratosMPICore/AssembleTestSuites.md)
- [AssembleTestSuites](../../../../functions/tests/test_KratosCore/AssembleTestSuites.md)