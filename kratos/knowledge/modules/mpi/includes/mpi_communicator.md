---
type: C++ Module
title: mpi_communicator
resource: mpi/includes/mpi_communicator.h#L1-L2239
generated:
  by: okf-rs/0.7.0
relationships:
  imports:
  - target: external/string
    resolved_by: tree-sitter
    confidence: exact
  - target: external/iostream
    resolved_by: tree-sitter
    confidence: exact
  - target: external/sstream
    resolved_by: tree-sitter
    confidence: exact
  - target: external/cstddef
    resolved_by: tree-sitter
    confidence: exact
  - target: external/mpi-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/includes-data-communicator-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/includes-define-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/includes-model-part-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/includes-mpi-serializer-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/mpi-includes-mpi-data-communicator-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/utilities-timer-h
    resolved_by: tree-sitter
    confidence: exact
  member_of:
  - target: packages/KratosMultiphysics
    resolved_by: tree-sitter
    confidence: exact
---

# Contains

- [GetMessageSize](../../../functions/mpi/includes/mpi_communicator/SendTraits< Vector >/GetMessageSize.md)
- [GetMessageSize](../../../functions/mpi/includes/mpi_communicator/SendTraits< Matrix >/GetMessageSize.md)
- [GetMessageSize](../../../functions/mpi/includes/mpi_communicator/SendTraits< DenseVector<TVectorValue> >/GetMessageSize.md)
- [GetMessageSize](../../../functions/mpi/includes/mpi_communicator/SendTraits< Kratos/VariablesListDataValueContainer >/GetMessageSize.md)
- [GetMessageSize](../../../functions/mpi/includes/mpi_communicator/SendTraits< Node/DofsContainerType >/GetMessageSize.md)
- [DirectCopyTransfer](../../../classes/mpi/includes/mpi_communicator/DirectCopyTransfer.md)
- [WriteBuffer](../../../functions/mpi/includes/mpi_communicator/DirectCopyTransfer/WriteBuffer.md)
- [ReadBuffer](../../../functions/mpi/includes/mpi_communicator/DirectCopyTransfer/ReadBuffer.md)
- [DynamicArrayTypeTransfer](../../../classes/mpi/includes/mpi_communicator/DynamicArrayTypeTransfer.md)
- [WriteBuffer](../../../functions/mpi/includes/mpi_communicator/DynamicArrayTypeTransfer/WriteBuffer.md)
- [ReadBuffer](../../../functions/mpi/includes/mpi_communicator/DynamicArrayTypeTransfer/ReadBuffer.md)
- [WriteBuffer](../../../functions/mpi/includes/mpi_communicator/SendTools< Quaternion<double> >/WriteBuffer.md)
- [ReadBuffer](../../../functions/mpi/includes/mpi_communicator/SendTools< Quaternion<double> >/ReadBuffer.md)
- [WriteBuffer](../../../functions/mpi/includes/mpi_communicator/SendTools< DenseVector<TVectorValue> >/WriteBuffer.md)
- [ReadBuffer](../../../functions/mpi/includes/mpi_communicator/SendTools< DenseVector<TVectorValue> >/ReadBuffer.md)
- [WriteBuffer](../../../functions/mpi/includes/mpi_communicator/SendTools< Kratos/VariablesListDataValueContainer >/WriteBuffer.md)
- [ReadBuffer](../../../functions/mpi/includes/mpi_communicator/SendTools< Kratos/VariablesListDataValueContainer >/ReadBuffer.md)
- [WriteBuffer](../../../functions/mpi/includes/mpi_communicator/SendTools< Node/DofsContainerType >/WriteBuffer.md)
- [ReadBuffer](../../../functions/mpi/includes/mpi_communicator/SendTools< Node/DofsContainerType >/ReadBuffer.md)
- [BufferAllocation](../../../classes/mpi/includes/mpi_communicator/BufferAllocation.md)
- [GetSendSize](../../../functions/mpi/includes/mpi_communicator/BufferAllocation<TDatabaseAccess, true>/GetSendSize.md)
- [GetSendSize](../../../functions/mpi/includes/mpi_communicator/BufferAllocation<TDatabaseAccess, true>/GetSendSize-2.md)
- [GetSendSize](../../../functions/mpi/includes/mpi_communicator/BufferAllocation<TDatabaseAccess, false>/GetSendSize.md)
- [GetSendSize](../../../functions/mpi/includes/mpi_communicator/BufferAllocation<TDatabaseAccess, false>/GetSendSize-2.md)
- [NodalContainerAccess](../../../classes/mpi/includes/mpi_communicator/NodalContainerAccess.md)
- [ElementalContainerAccess](../../../classes/mpi/includes/mpi_communicator/ElementalContainerAccess.md)
- [NodalSolutionStepValueAccess](../../../classes/mpi/includes/mpi_communicator/NodalSolutionStepValueAccess.md)
- [NodalSolutionStepValueAccess](../../../functions/mpi/includes/mpi_communicator/NodalSolutionStepValueAccess/NodalSolutionStepValueAccess.md)
- [NodalDataAccess](../../../classes/mpi/includes/mpi_communicator/NodalDataAccess.md)
- [NodalDataAccess](../../../functions/mpi/includes/mpi_communicator/NodalDataAccess/NodalDataAccess.md)
- [NodalFlagsAccess](../../../classes/mpi/includes/mpi_communicator/NodalFlagsAccess.md)
- [NodalFlagsAccess](../../../functions/mpi/includes/mpi_communicator/NodalFlagsAccess/NodalFlagsAccess.md)
- [NodalSolutionStepDataAccess](../../../classes/mpi/includes/mpi_communicator/NodalSolutionStepDataAccess.md)
- [DofIdAccess](../../../classes/mpi/includes/mpi_communicator/DofIdAccess.md)
- [ElementalDataAccess](../../../classes/mpi/includes/mpi_communicator/ElementalDataAccess.md)
- [ElementalDataAccess](../../../functions/mpi/includes/mpi_communicator/ElementalDataAccess/ElementalDataAccess.md)
- [ElementalFlagsAccess](../../../classes/mpi/includes/mpi_communicator/ElementalFlagsAccess.md)
- [ElementalFlagsAccess](../../../functions/mpi/includes/mpi_communicator/ElementalFlagsAccess/ElementalFlagsAccess.md)
- [MeshAccess](../../../classes/mpi/includes/mpi_communicator/MeshAccess.md)
- [Operation](../../../classes/mpi/includes/mpi_communicator/Operation.md)
- [MPICommunicator](../../../functions/mpi/includes/mpi_communicator/MPICommunicator.md)
- [MPICommunicator](../../../functions/mpi/includes/mpi_communicator/MPICommunicator-2.md)
- [MPICommunicator](../../../functions/mpi/includes/mpi_communicator/MPICommunicator-3.md)
- [Create](../../../functions/mpi/includes/mpi_communicator/Create.md)
- [IsDistributed](../../../functions/mpi/includes/mpi_communicator/IsDistributed.md)
- [SynchronizeNodalSolutionStepsData](../../../functions/mpi/includes/mpi_communicator/SynchronizeNodalSolutionStepsData.md)
- [SynchronizeDofs](../../../functions/mpi/includes/mpi_communicator/SynchronizeDofs.md)
- [SynchronizeVariable](../../../functions/mpi/includes/mpi_communicator/SynchronizeVariable.md)
- [SynchronizeVariable](../../../functions/mpi/includes/mpi_communicator/SynchronizeVariable-2.md)
- [SynchronizeVariable](../../../functions/mpi/includes/mpi_communicator/SynchronizeVariable-3.md)
- [SynchronizeVariable](../../../functions/mpi/includes/mpi_communicator/SynchronizeVariable-4.md)
- [SynchronizeVariable](../../../functions/mpi/includes/mpi_communicator/SynchronizeVariable-5.md)
- [SynchronizeVariable](../../../functions/mpi/includes/mpi_communicator/SynchronizeVariable-6.md)
- [SynchronizeVariable](../../../functions/mpi/includes/mpi_communicator/SynchronizeVariable-7.md)
- [SynchronizeVariable](../../../functions/mpi/includes/mpi_communicator/SynchronizeVariable-8.md)
- [SynchronizeVariable](../../../functions/mpi/includes/mpi_communicator/SynchronizeVariable-9.md)
- [SynchronizeVariable](../../../functions/mpi/includes/mpi_communicator/SynchronizeVariable-10.md)
- [SynchronizeNonHistoricalVariable](../../../functions/mpi/includes/mpi_communicator/SynchronizeNonHistoricalVariable.md)
- [SynchronizeNonHistoricalVariable](../../../functions/mpi/includes/mpi_communicator/SynchronizeNonHistoricalVariable-2.md)
- [SynchronizeNonHistoricalVariable](../../../functions/mpi/includes/mpi_communicator/SynchronizeNonHistoricalVariable-3.md)
- [SynchronizeNonHistoricalVariable](../../../functions/mpi/includes/mpi_communicator/SynchronizeNonHistoricalVariable-4.md)
- [SynchronizeNonHistoricalVariable](../../../functions/mpi/includes/mpi_communicator/SynchronizeNonHistoricalVariable-5.md)
- [SynchronizeNonHistoricalVariable](../../../functions/mpi/includes/mpi_communicator/SynchronizeNonHistoricalVariable-6.md)
- [SynchronizeNonHistoricalVariable](../../../functions/mpi/includes/mpi_communicator/SynchronizeNonHistoricalVariable-7.md)
- [SynchronizeNonHistoricalVariable](../../../functions/mpi/includes/mpi_communicator/SynchronizeNonHistoricalVariable-8.md)
- [SynchronizeNonHistoricalVariable](../../../functions/mpi/includes/mpi_communicator/SynchronizeNonHistoricalVariable-9.md)
- [SynchronizeNonHistoricalVariable](../../../functions/mpi/includes/mpi_communicator/SynchronizeNonHistoricalVariable-10.md)
- [SynchronizeCurrentDataToMax](../../../functions/mpi/includes/mpi_communicator/SynchronizeCurrentDataToMax.md)
- [SynchronizeNonHistoricalDataToMax](../../../functions/mpi/includes/mpi_communicator/SynchronizeNonHistoricalDataToMax.md)
- [SynchronizeCurrentDataToAbsMax](../../../functions/mpi/includes/mpi_communicator/SynchronizeCurrentDataToAbsMax.md)
- [SynchronizeNonHistoricalDataToAbsMax](../../../functions/mpi/includes/mpi_communicator/SynchronizeNonHistoricalDataToAbsMax.md)
- [SynchronizeCurrentDataToMin](../../../functions/mpi/includes/mpi_communicator/SynchronizeCurrentDataToMin.md)
- [SynchronizeNonHistoricalDataToMin](../../../functions/mpi/includes/mpi_communicator/SynchronizeNonHistoricalDataToMin.md)
- [SynchronizeCurrentDataToAbsMin](../../../functions/mpi/includes/mpi_communicator/SynchronizeCurrentDataToAbsMin.md)
- [SynchronizeNonHistoricalDataToAbsMin](../../../functions/mpi/includes/mpi_communicator/SynchronizeNonHistoricalDataToAbsMin.md)
- [AssembleCurrentData](../../../functions/mpi/includes/mpi_communicator/AssembleCurrentData.md)
- [AssembleCurrentData](../../../functions/mpi/includes/mpi_communicator/AssembleCurrentData-2.md)
- [AssembleCurrentData](../../../functions/mpi/includes/mpi_communicator/AssembleCurrentData-3.md)
- [AssembleCurrentData](../../../functions/mpi/includes/mpi_communicator/AssembleCurrentData-4.md)
- [AssembleCurrentData](../../../functions/mpi/includes/mpi_communicator/AssembleCurrentData-5.md)
- [AssembleNonHistoricalData](../../../functions/mpi/includes/mpi_communicator/AssembleNonHistoricalData.md)
- [AssembleNonHistoricalData](../../../functions/mpi/includes/mpi_communicator/AssembleNonHistoricalData-2.md)
- [AssembleNonHistoricalData](../../../functions/mpi/includes/mpi_communicator/AssembleNonHistoricalData-3.md)
- [AssembleNonHistoricalData](../../../functions/mpi/includes/mpi_communicator/AssembleNonHistoricalData-4.md)
- [AssembleNonHistoricalData](../../../functions/mpi/includes/mpi_communicator/AssembleNonHistoricalData-5.md)
- [AssembleNonHistoricalData](../../../functions/mpi/includes/mpi_communicator/AssembleNonHistoricalData-6.md)
- [SynchronizeElementalNonHistoricalVariable](../../../functions/mpi/includes/mpi_communicator/SynchronizeElementalNonHistoricalVariable.md)
- [SynchronizeElementalNonHistoricalVariable](../../../functions/mpi/includes/mpi_communicator/SynchronizeElementalNonHistoricalVariable-2.md)
- [SynchronizeElementalNonHistoricalVariable](../../../functions/mpi/includes/mpi_communicator/SynchronizeElementalNonHistoricalVariable-3.md)
- [SynchronizeElementalNonHistoricalVariable](../../../functions/mpi/includes/mpi_communicator/SynchronizeElementalNonHistoricalVariable-4.md)
- [SynchronizeElementalNonHistoricalVariable](../../../functions/mpi/includes/mpi_communicator/SynchronizeElementalNonHistoricalVariable-5.md)
- [SynchronizeElementalNonHistoricalVariable](../../../functions/mpi/includes/mpi_communicator/SynchronizeElementalNonHistoricalVariable-6.md)
- [SynchronizeElementalNonHistoricalVariable](../../../functions/mpi/includes/mpi_communicator/SynchronizeElementalNonHistoricalVariable-7.md)
- [TransferObjects](../../../functions/mpi/includes/mpi_communicator/TransferObjects.md)
- [TransferObjects](../../../functions/mpi/includes/mpi_communicator/TransferObjects-2.md)
- [TransferObjects](../../../functions/mpi/includes/mpi_communicator/TransferObjects-3.md)
- [TransferObjects](../../../functions/mpi/includes/mpi_communicator/TransferObjects-4.md)
- [TransferObjects](../../../functions/mpi/includes/mpi_communicator/TransferObjects-5.md)
- [TransferObjects](../../../functions/mpi/includes/mpi_communicator/TransferObjects-6.md)
- [SynchronizeOrNodalFlags](../../../functions/mpi/includes/mpi_communicator/SynchronizeOrNodalFlags.md)
- [SynchronizeAndNodalFlags](../../../functions/mpi/includes/mpi_communicator/SynchronizeAndNodalFlags.md)
- [SynchronizeNodalFlags](../../../functions/mpi/includes/mpi_communicator/SynchronizeNodalFlags.md)
- [SynchronizeElementalFlags](../../../functions/mpi/includes/mpi_communicator/SynchronizeElementalFlags.md)
- [Info](../../../functions/mpi/includes/mpi_communicator/Info.md)
- [PrintNodesId](../../../functions/mpi/includes/mpi_communicator/PrintNodesId.md)
- [PrintNodesId](../../../functions/mpi/includes/mpi_communicator/PrintNodesId-2.md)
- [AsyncSendAndReceiveObjects](../../../functions/mpi/includes/mpi_communicator/AsyncSendAndReceiveObjects.md)
- [SynchronizeFixedSizeValues](../../../functions/mpi/includes/mpi_communicator/SynchronizeFixedSizeValues.md)
- [SynchronizeDynamicVectorValues](../../../functions/mpi/includes/mpi_communicator/SynchronizeDynamicVectorValues.md)
- [SynchronizeDynamicMatrixValues](../../../functions/mpi/includes/mpi_communicator/SynchronizeDynamicMatrixValues.md)
- [AssembleFixedSizeValues](../../../functions/mpi/includes/mpi_communicator/AssembleFixedSizeValues.md)
- [AssembleDynamicVectorValues](../../../functions/mpi/includes/mpi_communicator/AssembleDynamicVectorValues.md)
- [AssembleDynamicMatrixValues](../../../functions/mpi/includes/mpi_communicator/AssembleDynamicMatrixValues.md)
- [ReduceValues](../../../functions/mpi/includes/mpi_communicator/ReduceValues.md)
- [ReduceValues](../../../functions/mpi/includes/mpi_communicator/ReduceValues-2.md)
- [ReduceValues](../../../functions/mpi/includes/mpi_communicator/ReduceValues-3.md)
- [ReduceValues](../../../functions/mpi/includes/mpi_communicator/ReduceValues-4.md)
- [ReduceValues](../../../functions/mpi/includes/mpi_communicator/ReduceValues-5.md)
- [ReduceValues](../../../functions/mpi/includes/mpi_communicator/ReduceValues-6.md)
- [ReduceValues](../../../functions/mpi/includes/mpi_communicator/ReduceValues-7.md)
- [ReduceValues](../../../functions/mpi/includes/mpi_communicator/ReduceValues-8.md)
- [ReduceValues](../../../functions/mpi/includes/mpi_communicator/ReduceValues-9.md)
- [TransferDistributedValues](../../../functions/mpi/includes/mpi_communicator/TransferDistributedValues.md)
- [TransferDistributedValuesUnknownSize](../../../functions/mpi/includes/mpi_communicator/TransferDistributedValuesUnknownSize.md)
- [AllocateBuffer](../../../functions/mpi/includes/mpi_communicator/AllocateBuffer.md)
- [FillBuffer](../../../functions/mpi/includes/mpi_communicator/FillBuffer.md)
- [FillBuffer](../../../functions/mpi/includes/mpi_communicator/FillBuffer-2.md)
- [UpdateValues](../../../functions/mpi/includes/mpi_communicator/UpdateValues.md)
- [UpdateValues](../../../functions/mpi/includes/mpi_communicator/UpdateValues-2.md)
- [MatchDynamicVectorSizes](../../../functions/mpi/includes/mpi_communicator/MatchDynamicVectorSizes.md)
- [MatchDynamicMatrixSizes](../../../functions/mpi/includes/mpi_communicator/MatchDynamicMatrixSizes.md)

# Imports

- `<string>`
- `<iostream>`
- `<sstream>`
- `<cstddef>`
- `"mpi.h"`
- `"includes/data_communicator.h"`
- `"includes/define.h"`
- `"includes/model_part.h"`
- `"includes/mpi_serializer.h"`
- `"mpi/includes/mpi_data_communicator.h"`
- `"utilities/timer.h"`

# Member of

- [KratosMultiphysics](../../../packages/KratosMultiphysics.md)