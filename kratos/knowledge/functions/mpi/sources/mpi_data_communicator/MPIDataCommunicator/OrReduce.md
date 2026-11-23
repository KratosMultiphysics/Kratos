---
type: C++ Method
title: OrReduce
resource: mpi/sources/mpi_data_communicator.cpp#L364-L367
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

`bool MPIDataCommunicator::OrReduce(const bool Value, const int Root) const`

# Calls

- [GetDefined](../../../../../functions/containers/flags/Flags/GetDefined.md)
- [GetFlags](../../../../../functions/containers/flags/Flags/GetFlags.md)
- [SetDefined](../../../../../functions/containers/flags/Flags/SetDefined.md)