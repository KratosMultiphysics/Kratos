---
type: C++ Function
title: MatchDynamicVectorSizes
resource: mpi/includes/mpi_communicator.h#L2010-L2095
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/mpi/includes/mpi_communicator/SynchronizeDynamicVectorValues
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/mpi/includes/mpi_communicator/AssembleDynamicVectorValues
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void MatchDynamicVectorSizes( TSourceAccess SourceType, TDestinationAccess DestinationType, TDatabaseAccess& rAccess)`

# Called by

- [SynchronizeDynamicVectorValues](../../../../functions/mpi/includes/mpi_communicator/SynchronizeDynamicVectorValues.md)
- [AssembleDynamicVectorValues](../../../../functions/mpi/includes/mpi_communicator/AssembleDynamicVectorValues.md)