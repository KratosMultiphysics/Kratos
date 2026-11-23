---
type: C++ Method
title: AndReduceAll
resource: mpi/sources/mpi_data_communicator.cpp#L387-L390
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/containers/flags/Flags/GetDefined
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/containers/flags/Flags/GetFlags
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/containers/flags/Flags/SetDefined
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`bool MPIDataCommunicator::AndReduceAll(const bool Value) const`

# Calls

- [GetDefined](../../../../../functions/containers/flags/Flags/GetDefined.md)
- [GetFlags](../../../../../functions/containers/flags/Flags/GetFlags.md)
- [SetDefined](../../../../../functions/containers/flags/Flags/SetDefined.md)