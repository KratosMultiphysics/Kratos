---
type: C++ Function
title: TotalSize
resource: containers/variables_list_data_value_container.h#L392-L398
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/containers/variables_list_data_value_container/CloneFront
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/containers/variables_list_data_value_container/PushFront
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/containers/variables_list_data_value_container/Position
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/mpi/includes/mpi_communicator/SendTraits< Kratos/VariablesListDataValueContainer >/GetMessageSize
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/mpi/includes/mpi_communicator/SendTools< Kratos/VariablesListDataValueContainer >/WriteBuffer
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/mpi/includes/mpi_communicator/SendTools< Kratos/VariablesListDataValueContainer >/ReadBuffer
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`SizeType TotalSize() const`

# Called by

- [CloneFront](../../../functions/containers/variables_list_data_value_container/CloneFront.md)
- [PushFront](../../../functions/containers/variables_list_data_value_container/PushFront.md)
- [Position](../../../functions/containers/variables_list_data_value_container/Position.md)
- [GetMessageSize](../../../functions/mpi/includes/mpi_communicator/SendTraits< Kratos/VariablesListDataValueContainer >/GetMessageSize.md)
- [WriteBuffer](../../../functions/mpi/includes/mpi_communicator/SendTools< Kratos/VariablesListDataValueContainer >/WriteBuffer.md)
- [ReadBuffer](../../../functions/mpi/includes/mpi_communicator/SendTools< Kratos/VariablesListDataValueContainer >/ReadBuffer.md)