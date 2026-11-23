---
type: Python Method
title: testGatherOnMaster
resource: mpi/tests/test_gather_modelpart_utility.py#L106-L119
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/mpi/tests/test_gather_modelpart_utility/SetPressureOnLocalNodes
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/KratosUnittest/isclose
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`def testGatherOnMaster(self):`

# Calls

- [SetPressureOnLocalNodes](../../../../../functions/mpi/tests/test_gather_modelpart_utility/SetPressureOnLocalNodes.md)
- [isclose](../../../../../functions/python_scripts/KratosUnittest/isclose.md)