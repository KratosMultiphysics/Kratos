---
type: C++ Method
title: SetUpMPIEnvironmentDetail
resource: sources/parallel_environment.cpp#L274-L280
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/sources/parallel_environment/ParallelEnvironment/MPIIsInitialized
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/parallel_environment/ParallelEnvironment/MPIIsFinalized
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/sources/parallel_environment/ParallelEnvironment/SetUpMPIEnvironment
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void ParallelEnvironment::SetUpMPIEnvironmentDetail(EnvironmentManager::Pointer pEnvironmentManager)`

# Calls

- [MPIIsInitialized](../../../../functions/sources/parallel_environment/ParallelEnvironment/MPIIsInitialized.md)
- [MPIIsFinalized](../../../../functions/sources/parallel_environment/ParallelEnvironment/MPIIsFinalized.md)

# Called by

- [SetUpMPIEnvironment](../../../../functions/sources/parallel_environment/ParallelEnvironment/SetUpMPIEnvironment.md)