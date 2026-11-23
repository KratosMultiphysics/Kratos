---
type: C++ Function
title: InitializeMPIParallelRun
resource: mpi/python/kratos_mpi_python.cpp#L32-L62
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
  called_by:
  - target: functions/python_interface/__ModuleInitDetail
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void InitializeMPIParallelRun()`

# Calls

- [RegisterDataCommunicator](../../../../functions/sources/parallel_environment/ParallelEnvironment/RegisterDataCommunicator.md)
- [HasDataCommunicator](../../../../functions/sources/parallel_environment/ParallelEnvironment/HasDataCommunicator.md)

# Called by

- [__ModuleInitDetail](../../../../functions/python_interface/__ModuleInitDetail.md)