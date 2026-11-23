---
type: C++ Method
title: MPIIsFinalized
resource: sources/parallel_environment.cpp#L147-L151
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/sources/parallel_environment/ParallelEnvironment/MPIIsFinalizedDetail
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/sources/parallel_environment/ParallelEnvironment/SetUpMPIEnvironmentDetail
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`bool ParallelEnvironment::MPIIsFinalized()`

# Calls

- [MPIIsFinalizedDetail](../../../../functions/sources/parallel_environment/ParallelEnvironment/MPIIsFinalizedDetail.md)

# Called by

- [SetUpMPIEnvironmentDetail](../../../../functions/sources/parallel_environment/ParallelEnvironment/SetUpMPIEnvironmentDetail.md)