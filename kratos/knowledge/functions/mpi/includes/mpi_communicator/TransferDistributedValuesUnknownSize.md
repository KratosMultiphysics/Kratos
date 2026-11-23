---
type: C++ Function
title: TransferDistributedValuesUnknownSize
resource: mpi/includes/mpi_communicator.h#L1866-L1911
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/mpi/includes/mpi_communicator/SynchronizeNodalSolutionStepsData
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void TransferDistributedValuesUnknownSize( TSourceAccess SourceType, TDestinationAccess DestinationType, TDatabaseAccess& rAccess, TReductionOperation Reduction)`

# Called by

- [SynchronizeNodalSolutionStepsData](../../../../functions/mpi/includes/mpi_communicator/SynchronizeNodalSolutionStepsData.md)