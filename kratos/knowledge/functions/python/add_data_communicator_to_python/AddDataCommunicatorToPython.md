---
type: C++ Function
title: AddDataCommunicatorToPython
resource: python/add_data_communicator_to_python.cpp#L136-L175
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/includes/data_communicator/Broadcast
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/python/kratos_python/PYBIND11_MODULE
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void AddDataCommunicatorToPython(pybind11::module &m)`

# Calls

- [Broadcast](../../../functions/includes/data_communicator/Broadcast.md)

# Called by

- [PYBIND11_MODULE](../../../functions/python/kratos_python/PYBIND11_MODULE.md)