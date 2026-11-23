---
type: Python Method
title: testScatterFromMaster
resource: mpi/tests/test_gather_modelpart_utility.py#L121-L134
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/mpi/tests/test_gather_modelpart_utility/SetPressureOnAllNodes
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/mpi/utilities/gather_modelpart_utility/GatherModelPartUtility/ScatterFromMaster
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/KratosUnittest/isclose
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`def testScatterFromMaster(self):`

# Calls

- [SetPressureOnAllNodes](../../../../../functions/mpi/tests/test_gather_modelpart_utility/SetPressureOnAllNodes.md)
- [ScatterFromMaster](../../../../../functions/mpi/utilities/gather_modelpart_utility/GatherModelPartUtility/ScatterFromMaster.md)
- [isclose](../../../../../functions/python_scripts/KratosUnittest/isclose.md)