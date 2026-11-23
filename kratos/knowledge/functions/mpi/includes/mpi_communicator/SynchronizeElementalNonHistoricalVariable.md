---
type: C++ Function
title: SynchronizeElementalNonHistoricalVariable
resource: mpi/includes/mpi_communicator.h#L1141-L1146
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/mpi/includes/mpi_communicator/SynchronizeFixedSizeValues
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/mpi/includes/mpi_communicator/AssembleDynamicVectorValues
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/mpi/includes/mpi_communicator/SynchronizeDynamicVectorValues
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/mpi/includes/mpi_communicator/SynchronizeDynamicMatrixValues
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`bool SynchronizeElementalNonHistoricalVariable(Variable<int> const& ThisVariable) override`

# Calls

- [SynchronizeFixedSizeValues](../../../../functions/mpi/includes/mpi_communicator/SynchronizeFixedSizeValues.md)
- [AssembleDynamicVectorValues](../../../../functions/mpi/includes/mpi_communicator/AssembleDynamicVectorValues.md)
- [SynchronizeDynamicVectorValues](../../../../functions/mpi/includes/mpi_communicator/SynchronizeDynamicVectorValues.md)
- [SynchronizeDynamicMatrixValues](../../../../functions/mpi/includes/mpi_communicator/SynchronizeDynamicMatrixValues.md)