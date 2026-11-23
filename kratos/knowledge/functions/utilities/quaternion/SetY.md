---
type: C++ Function
title: SetY
resource: utilities/quaternion.h#L130
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/mpi/includes/mpi_communicator/SendTools< Quaternion<double> >/ReadBuffer
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python/add_quaternion_to_python/QuaternionSetY
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/quaternion/Quaternion
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`inline void SetY(const T& value)`

# Called by

- [ReadBuffer](../../../functions/mpi/includes/mpi_communicator/SendTools< Quaternion<double> >/ReadBuffer.md)
- [QuaternionSetY](../../../functions/python/add_quaternion_to_python/QuaternionSetY.md)
- [Quaternion](../../../functions/utilities/quaternion/Quaternion.md)