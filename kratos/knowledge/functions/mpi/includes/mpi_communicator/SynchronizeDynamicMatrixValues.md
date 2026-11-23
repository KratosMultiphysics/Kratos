---
type: C++ Function
title: SynchronizeDynamicMatrixValues
resource: mpi/includes/mpi_communicator.h#L1587-L1598
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/mpi/includes/mpi_communicator/MatchDynamicMatrixSizes
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

`void SynchronizeDynamicMatrixValues(TDatabaseAccess& rVariableAccess)`

# Calls

- [MatchDynamicMatrixSizes](../../../../functions/mpi/includes/mpi_communicator/MatchDynamicMatrixSizes.md)
- [TransferDistributedValues](../../../../functions/mpi/includes/mpi_communicator/TransferDistributedValues.md)

# Called by

- [SynchronizeVariable](../../../../functions/mpi/includes/mpi_communicator/SynchronizeVariable.md)
- [SynchronizeNonHistoricalVariable](../../../../functions/mpi/includes/mpi_communicator/SynchronizeNonHistoricalVariable.md)
- [SynchronizeElementalNonHistoricalVariable](../../../../functions/mpi/includes/mpi_communicator/SynchronizeElementalNonHistoricalVariable.md)