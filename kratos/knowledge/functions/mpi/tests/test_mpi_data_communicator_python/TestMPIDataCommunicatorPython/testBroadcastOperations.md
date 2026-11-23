---
type: Python Method
title: testBroadcastOperations
resource: mpi/tests/test_mpi_data_communicator_python.py#L177-L195
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/includes/data_communicator/Broadcast
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/KratosUnittest/TestCase/assertVectorAlmostEqual
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/KratosUnittest/TestCase/assertMatrixAlmostEqual
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`def testBroadcastOperations(self):`

# Calls

- [Broadcast](../../../../../functions/includes/data_communicator/Broadcast.md)
- [assertVectorAlmostEqual](../../../../../functions/python_scripts/KratosUnittest/TestCase/assertVectorAlmostEqual.md)
- [assertMatrixAlmostEqual](../../../../../functions/python_scripts/KratosUnittest/TestCase/assertMatrixAlmostEqual.md)