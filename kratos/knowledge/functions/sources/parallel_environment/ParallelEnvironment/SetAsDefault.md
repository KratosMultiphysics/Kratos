---
type: C++ Method
title: SetAsDefault
resource: sources/parallel_environment.cpp#L259-L266
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/sources/parallel_environment/ParallelEnvironment/RegisterDataCommunicatorDetail
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/parallel_environment/ParallelEnvironment/SetDefaultDataCommunicatorDetail
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void ParallelEnvironment::SetAsDefault( std::unordered_map<std::string, DataCommunicator::UniquePointer>::iterator& rThisCommunicator)`

# Called by

- [RegisterDataCommunicatorDetail](../../../../functions/sources/parallel_environment/ParallelEnvironment/RegisterDataCommunicatorDetail.md)
- [SetDefaultDataCommunicatorDetail](../../../../functions/sources/parallel_environment/ParallelEnvironment/SetDefaultDataCommunicatorDetail.md)