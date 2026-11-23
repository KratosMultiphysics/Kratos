---
type: C++ Function
title: MatchDynamicMatrixSizes
resource: mpi/includes/mpi_communicator.h#L2101-L2188
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/mpi/includes/mpi_communicator/SynchronizeDynamicMatrixValues
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/mpi/includes/mpi_communicator/AssembleDynamicMatrixValues
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void MatchDynamicMatrixSizes( TSourceAccess SourceType, TDestinationAccess DestinationType, TDatabaseAccess& rAccess)`

# Called by

- [SynchronizeDynamicMatrixValues](../../../../functions/mpi/includes/mpi_communicator/SynchronizeDynamicMatrixValues.md)
- [AssembleDynamicMatrixValues](../../../../functions/mpi/includes/mpi_communicator/AssembleDynamicMatrixValues.md)