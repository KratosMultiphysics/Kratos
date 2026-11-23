---
type: C++ Function
title: AssembleDynamicVectorValues
resource: mpi/includes/mpi_communicator.h#L1616-L1636
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
  - target: functions/mpi/includes/mpi_communicator/AssembleCurrentData
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/mpi/includes/mpi_communicator/AssembleNonHistoricalData
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/mpi/includes/mpi_communicator/SynchronizeElementalNonHistoricalVariable
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void AssembleDynamicVectorValues(TDatabaseAccess& rVariableAccess)`

# Calls

- [MatchDynamicVectorSizes](../../../../functions/mpi/includes/mpi_communicator/MatchDynamicVectorSizes.md)
- [TransferDistributedValues](../../../../functions/mpi/includes/mpi_communicator/TransferDistributedValues.md)

# Called by

- [AssembleCurrentData](../../../../functions/mpi/includes/mpi_communicator/AssembleCurrentData.md)
- [AssembleNonHistoricalData](../../../../functions/mpi/includes/mpi_communicator/AssembleNonHistoricalData.md)
- [SynchronizeElementalNonHistoricalVariable](../../../../functions/mpi/includes/mpi_communicator/SynchronizeElementalNonHistoricalVariable.md)