---
type: C++ Method
title: RegisterDataCommunicatorDetail
resource: sources/parallel_environment.cpp#L282-L307
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/sources/parallel_environment/ParallelEnvironment/SetAsDefault
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/sources/parallel_environment/ParallelEnvironment/RegisterDataCommunicator
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/parallel_environment/ParallelEnvironment/ParallelEnvironment
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void ParallelEnvironment::RegisterDataCommunicatorDetail( const std::string& Name, DataCommunicator::UniquePointer pPrototype, const bool Default)`

# Calls

- [SetAsDefault](../../../../functions/sources/parallel_environment/ParallelEnvironment/SetAsDefault.md)

# Called by

- [RegisterDataCommunicator](../../../../functions/sources/parallel_environment/ParallelEnvironment/RegisterDataCommunicator.md)
- [ParallelEnvironment](../../../../functions/sources/parallel_environment/ParallelEnvironment/ParallelEnvironment.md)