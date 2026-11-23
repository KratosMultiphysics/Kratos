---
type: C++ Method
title: TearDown
resource: mpi/tests/test_utilities/mpi_test_environment.cpp#L65-L72
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/sources/parallel_environment/ParallelEnvironment/RegisterDataCommunicator
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/parallel_environment/ParallelEnvironment/UnregisterDataCommunicator
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void KratosMpiTestEnv::TearDown()`

# Calls

- [RegisterDataCommunicator](../../../../../../functions/sources/parallel_environment/ParallelEnvironment/RegisterDataCommunicator.md)
- [UnregisterDataCommunicator](../../../../../../functions/sources/parallel_environment/ParallelEnvironment/UnregisterDataCommunicator.md)