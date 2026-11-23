---
type: C++ Function
title: TransferDistributedValues
resource: mpi/includes/mpi_communicator.h#L1820-L1859
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/mpi/includes/mpi_communicator/AllocateBuffer
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/mpi/includes/mpi_communicator/SynchronizeDofs
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/mpi/includes/mpi_communicator/SynchronizeCurrentDataToMax
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/mpi/includes/mpi_communicator/SynchronizeNonHistoricalDataToMax
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/mpi/includes/mpi_communicator/SynchronizeCurrentDataToAbsMax
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/mpi/includes/mpi_communicator/SynchronizeNonHistoricalDataToAbsMax
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/mpi/includes/mpi_communicator/SynchronizeCurrentDataToMin
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/mpi/includes/mpi_communicator/SynchronizeNonHistoricalDataToMin
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/mpi/includes/mpi_communicator/SynchronizeCurrentDataToAbsMin
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/mpi/includes/mpi_communicator/SynchronizeNonHistoricalDataToAbsMin
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/mpi/includes/mpi_communicator/SynchronizeOrNodalFlags
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/mpi/includes/mpi_communicator/SynchronizeAndNodalFlags
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/mpi/includes/mpi_communicator/SynchronizeNodalFlags
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/mpi/includes/mpi_communicator/SynchronizeElementalFlags
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/mpi/includes/mpi_communicator/SynchronizeFixedSizeValues
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/mpi/includes/mpi_communicator/SynchronizeDynamicVectorValues
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/mpi/includes/mpi_communicator/SynchronizeDynamicMatrixValues
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/mpi/includes/mpi_communicator/AssembleFixedSizeValues
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/mpi/includes/mpi_communicator/AssembleDynamicVectorValues
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/mpi/includes/mpi_communicator/AssembleDynamicMatrixValues
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void TransferDistributedValues( TSourceAccess SourceType, TDestinationAccess DestinationType, TDatabaseAccess& rAccess, TReductionOperation Reduction)`

# Calls

- [AllocateBuffer](../../../../functions/mpi/includes/mpi_communicator/AllocateBuffer.md)

# Called by

- [SynchronizeDofs](../../../../functions/mpi/includes/mpi_communicator/SynchronizeDofs.md)
- [SynchronizeCurrentDataToMax](../../../../functions/mpi/includes/mpi_communicator/SynchronizeCurrentDataToMax.md)
- [SynchronizeNonHistoricalDataToMax](../../../../functions/mpi/includes/mpi_communicator/SynchronizeNonHistoricalDataToMax.md)
- [SynchronizeCurrentDataToAbsMax](../../../../functions/mpi/includes/mpi_communicator/SynchronizeCurrentDataToAbsMax.md)
- [SynchronizeNonHistoricalDataToAbsMax](../../../../functions/mpi/includes/mpi_communicator/SynchronizeNonHistoricalDataToAbsMax.md)
- [SynchronizeCurrentDataToMin](../../../../functions/mpi/includes/mpi_communicator/SynchronizeCurrentDataToMin.md)
- [SynchronizeNonHistoricalDataToMin](../../../../functions/mpi/includes/mpi_communicator/SynchronizeNonHistoricalDataToMin.md)
- [SynchronizeCurrentDataToAbsMin](../../../../functions/mpi/includes/mpi_communicator/SynchronizeCurrentDataToAbsMin.md)
- [SynchronizeNonHistoricalDataToAbsMin](../../../../functions/mpi/includes/mpi_communicator/SynchronizeNonHistoricalDataToAbsMin.md)
- [SynchronizeOrNodalFlags](../../../../functions/mpi/includes/mpi_communicator/SynchronizeOrNodalFlags.md)
- [SynchronizeAndNodalFlags](../../../../functions/mpi/includes/mpi_communicator/SynchronizeAndNodalFlags.md)
- [SynchronizeNodalFlags](../../../../functions/mpi/includes/mpi_communicator/SynchronizeNodalFlags.md)
- [SynchronizeElementalFlags](../../../../functions/mpi/includes/mpi_communicator/SynchronizeElementalFlags.md)
- [SynchronizeFixedSizeValues](../../../../functions/mpi/includes/mpi_communicator/SynchronizeFixedSizeValues.md)
- [SynchronizeDynamicVectorValues](../../../../functions/mpi/includes/mpi_communicator/SynchronizeDynamicVectorValues.md)
- [SynchronizeDynamicMatrixValues](../../../../functions/mpi/includes/mpi_communicator/SynchronizeDynamicMatrixValues.md)
- [AssembleFixedSizeValues](../../../../functions/mpi/includes/mpi_communicator/AssembleFixedSizeValues.md)
- [AssembleDynamicVectorValues](../../../../functions/mpi/includes/mpi_communicator/AssembleDynamicVectorValues.md)
- [AssembleDynamicMatrixValues](../../../../functions/mpi/includes/mpi_communicator/AssembleDynamicMatrixValues.md)