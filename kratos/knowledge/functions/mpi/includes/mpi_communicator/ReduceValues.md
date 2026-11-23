---
type: C++ Function
title: ReduceValues
resource: mpi/includes/mpi_communicator.h#L1672-L1686
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/containers/flags/AssignFlags
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`std::size_t ReduceValues( const typename TDatabaseAccess::SendType* pBuffer, TDatabaseAccess& rAccess, typename TDatabaseAccess::IteratorType ContainerIterator, Operation<OperationType::Replace>)`

# Calls

- [AssignFlags](../../../../functions/containers/flags/AssignFlags.md)