---
type: Python Method
title: tearDown
resource: mpi/tests/test_data_communicator_factory.py#L15-L19
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/sources/parallel_environment/ParallelEnvironment/SetDefaultDataCommunicator
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/parallel_environment/ParallelEnvironment/UnregisterDataCommunicator
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`def tearDown(self):`

# Calls

- [SetDefaultDataCommunicator](../../../../../functions/sources/parallel_environment/ParallelEnvironment/SetDefaultDataCommunicator.md)
- [UnregisterDataCommunicator](../../../../../functions/sources/parallel_environment/ParallelEnvironment/UnregisterDataCommunicator.md)