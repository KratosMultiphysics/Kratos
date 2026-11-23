---
type: C++ Function
title: SynchronizeDynamicVectorValues
resource: mpi/includes/mpi_communicator.h#L1573-L1584
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/mpi/includes/mpi_communicator/MatchDynamicVectorSizes
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/mpi/includes/mpi_communicator/TransferDistributedValues
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/mpi/includes/mpi_communicator/SynchronizeVariable
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/mpi/includes/mpi_communicator/SynchronizeNonHistoricalVariable
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/mpi/includes/mpi_communicator/SynchronizeElementalNonHistoricalVariable
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void SynchronizeDynamicVectorValues(TDatabaseAccess& rVariableAccess)`

# Calls

- [MatchDynamicVectorSizes](../../../../functions/mpi/includes/mpi_communicator/MatchDynamicVectorSizes.md)
- [TransferDistributedValues](../../../../functions/mpi/includes/mpi_communicator/TransferDistributedValues.md)

# Called by

- [SynchronizeVariable](../../../../functions/mpi/includes/mpi_communicator/SynchronizeVariable.md)
- [SynchronizeNonHistoricalVariable](../../../../functions/mpi/includes/mpi_communicator/SynchronizeNonHistoricalVariable.md)
- [SynchronizeElementalNonHistoricalVariable](../../../../functions/mpi/includes/mpi_communicator/SynchronizeElementalNonHistoricalVariable.md)