---
type: C++ Method
title: RetrieveRegisteredName
resource: sources/parallel_environment.cpp#L235-L243
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/containers/distributed_csr_matrix/DistributedCsrMatrix/save
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`std::string ParallelEnvironment::RetrieveRegisteredName(const DataCommunicator& rComm)`

# Called by

- [save](../../../../functions/containers/distributed_csr_matrix/DistributedCsrMatrix/save.md)