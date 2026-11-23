---
type: C++ Method
title: SetUp
resource: mpi/tests/test_utilities/mpi_test_environment.cpp#L32-L63
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/sources/parallel_environment/ParallelEnvironment/RegisterDataCommunicator
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/parallel_environment/ParallelEnvironment/HasDataCommunicator
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void KratosMpiTestEnv::SetUp()`

# Calls

- [RegisterDataCommunicator](../../../../../../functions/sources/parallel_environment/ParallelEnvironment/RegisterDataCommunicator.md)
- [HasDataCommunicator](../../../../../../functions/sources/parallel_environment/ParallelEnvironment/HasDataCommunicator.md)