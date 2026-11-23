---
type: C++ Function
title: SynchronizeNonHistoricalVariable
resource: mpi/includes/mpi_communicator.h#L856-L861
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/mpi/includes/mpi_communicator/SynchronizeFixedSizeValues
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

`bool SynchronizeNonHistoricalVariable(Variable<int> const& rThisVariable) override`

# Calls

- [SynchronizeFixedSizeValues](../../../../functions/mpi/includes/mpi_communicator/SynchronizeFixedSizeValues.md)
- [SynchronizeDynamicVectorValues](../../../../functions/mpi/includes/mpi_communicator/SynchronizeDynamicVectorValues.md)
- [SynchronizeDynamicMatrixValues](../../../../functions/mpi/includes/mpi_communicator/SynchronizeDynamicMatrixValues.md)