---
type: C++ Method
title: SetDefaultDataCommunicatorDetail
resource: sources/parallel_environment.cpp#L339-L347
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/sources/parallel_environment/ParallelEnvironment/SetAsDefault
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/sources/parallel_environment/ParallelEnvironment/SetDefaultDataCommunicator
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void ParallelEnvironment::SetDefaultDataCommunicatorDetail(const std::string& rName)`

# Calls

- [SetAsDefault](../../../../functions/sources/parallel_environment/ParallelEnvironment/SetAsDefault.md)

# Called by

- [SetDefaultDataCommunicator](../../../../functions/sources/parallel_environment/ParallelEnvironment/SetDefaultDataCommunicator.md)