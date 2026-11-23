---
type: C++ Module
title: mpi_data_communicator
resource: mpi/sources/mpi_data_communicator.cpp#L1-L1657
generated:
  by: okf-rs/0.7.0
relationships:
  imports:
  - target: external/algorithm
    resolved_by: tree-sitter
    confidence: exact
  - target: external/includes-parallel-environment-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/mpi-includes-mpi-data-communicator-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/mpi-includes-mpi-manager-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/mpi-includes-mpi-message-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/mpi-utilities-data-communicator-factory-h
    resolved_by: tree-sitter
    confidence: exact
  member_of:
  - target: packages/KratosMultiphysics
    resolved_by: tree-sitter
    confidence: exact
---

# Contains

- [MPIDataCommunicator](../../../functions/mpi/sources/mpi_data_communicator/MPIDataCommunicator/MPIDataCommunicator.md)
- [Create](../../../functions/mpi/sources/mpi_data_communicator/MPIDataCommunicator/Create.md)
- [Barrier](../../../functions/mpi/sources/mpi_data_communicator/MPIDataCommunicator/Barrier.md)
- [AndReduce](../../../functions/mpi/sources/mpi_data_communicator/MPIDataCommunicator/AndReduce.md)
- [AndReduce](../../../functions/mpi/sources/mpi_data_communicator/MPIDataCommunicator/AndReduce-2.md)
- [OrReduce](../../../functions/mpi/sources/mpi_data_communicator/MPIDataCommunicator/OrReduce.md)
- [OrReduce](../../../functions/mpi/sources/mpi_data_communicator/MPIDataCommunicator/OrReduce-2.md)
- [AndReduceAll](../../../functions/mpi/sources/mpi_data_communicator/MPIDataCommunicator/AndReduceAll.md)
- [AndReduceAll](../../../functions/mpi/sources/mpi_data_communicator/MPIDataCommunicator/AndReduceAll-2.md)
- [OrReduceAll](../../../functions/mpi/sources/mpi_data_communicator/MPIDataCommunicator/OrReduceAll.md)
- [OrReduceAll](../../../functions/mpi/sources/mpi_data_communicator/MPIDataCommunicator/OrReduceAll-2.md)
- [GetMPICommunicator](../../../functions/mpi/sources/mpi_data_communicator/MPIDataCommunicator/GetMPICommunicator.md)
- [Rank](../../../functions/mpi/sources/mpi_data_communicator/MPIDataCommunicator/Rank.md)
- [Size](../../../functions/mpi/sources/mpi_data_communicator/MPIDataCommunicator/Size.md)
- [IsDistributed](../../../functions/mpi/sources/mpi_data_communicator/MPIDataCommunicator/IsDistributed.md)
- [IsDefinedOnThisRank](../../../functions/mpi/sources/mpi_data_communicator/MPIDataCommunicator/IsDefinedOnThisRank.md)
- [IsNullOnThisRank](../../../functions/mpi/sources/mpi_data_communicator/MPIDataCommunicator/IsNullOnThisRank.md)
- [Info](../../../functions/mpi/sources/mpi_data_communicator/MPIDataCommunicator/Info.md)
- [PrintInfo](../../../functions/mpi/sources/mpi_data_communicator/MPIDataCommunicator/PrintInfo.md)
- [PrintData](../../../functions/mpi/sources/mpi_data_communicator/MPIDataCommunicator/PrintData.md)
- [CheckMPIErrorCode](../../../functions/mpi/sources/mpi_data_communicator/MPIDataCommunicator/CheckMPIErrorCode.md)
- [BroadcastImpl](../../../functions/mpi/sources/mpi_data_communicator/MPIDataCommunicator/BroadcastImpl.md)
- [SendRecvImpl](../../../functions/mpi/sources/mpi_data_communicator/MPIDataCommunicator/SendRecvImpl.md)
- [SendRecvImpl](../../../functions/mpi/sources/mpi_data_communicator/MPIDataCommunicator/SendRecvImpl-2.md)
- [SendImpl](../../../functions/mpi/sources/mpi_data_communicator/MPIDataCommunicator/SendImpl.md)
- [RecvImpl](../../../functions/mpi/sources/mpi_data_communicator/MPIDataCommunicator/RecvImpl.md)
- [SynchronizeShapeDetail](../../../functions/mpi/sources/mpi_data_communicator/MPIDataCommunicator/SynchronizeShapeDetail.md)
- [SynchronizeShapeDetail](../../../functions/mpi/sources/mpi_data_communicator/MPIDataCommunicator/SynchronizeShapeDetail-2.md)
- [ReduceDetail](../../../functions/mpi/sources/mpi_data_communicator/MPIDataCommunicator/ReduceDetail.md)
- [ReduceDetail](../../../functions/mpi/sources/mpi_data_communicator/MPIDataCommunicator/ReduceDetail-2.md)
- [ReduceDetailVector](../../../functions/mpi/sources/mpi_data_communicator/MPIDataCommunicator/ReduceDetailVector.md)
- [AllReduceDetail](../../../functions/mpi/sources/mpi_data_communicator/MPIDataCommunicator/AllReduceDetail.md)
- [AllReduceDetail](../../../functions/mpi/sources/mpi_data_communicator/MPIDataCommunicator/AllReduceDetail-2.md)
- [AllReduceDetailVector](../../../functions/mpi/sources/mpi_data_communicator/MPIDataCommunicator/AllReduceDetailVector.md)
- [AllReduceDetailWithLocation](../../../functions/mpi/sources/mpi_data_communicator/MPIDataCommunicator/AllReduceDetailWithLocation.md)
- [ScanDetail](../../../functions/mpi/sources/mpi_data_communicator/MPIDataCommunicator/ScanDetail.md)
- [ScanDetail](../../../functions/mpi/sources/mpi_data_communicator/MPIDataCommunicator/ScanDetail-2.md)
- [ScanDetail](../../../functions/mpi/sources/mpi_data_communicator/MPIDataCommunicator/ScanDetail-3.md)
- [SendRecvDetail](../../../functions/mpi/sources/mpi_data_communicator/MPIDataCommunicator/SendRecvDetail.md)
- [SendRecvDetail](../../../functions/mpi/sources/mpi_data_communicator/MPIDataCommunicator/SendRecvDetail-2.md)
- [SendRecvDetail](../../../functions/mpi/sources/mpi_data_communicator/MPIDataCommunicator/SendRecvDetail-3.md)
- [SendDetail](../../../functions/mpi/sources/mpi_data_communicator/MPIDataCommunicator/SendDetail.md)
- [RecvDetail](../../../functions/mpi/sources/mpi_data_communicator/MPIDataCommunicator/RecvDetail.md)
- [BroadcastDetail](../../../functions/mpi/sources/mpi_data_communicator/MPIDataCommunicator/BroadcastDetail.md)
- [ScatterDetail](../../../functions/mpi/sources/mpi_data_communicator/MPIDataCommunicator/ScatterDetail.md)
- [ScatterDetail](../../../functions/mpi/sources/mpi_data_communicator/MPIDataCommunicator/ScatterDetail-2.md)
- [ScattervDetail](../../../functions/mpi/sources/mpi_data_communicator/MPIDataCommunicator/ScattervDetail.md)
- [ScattervDetail](../../../functions/mpi/sources/mpi_data_communicator/MPIDataCommunicator/ScattervDetail-2.md)
- [GatherDetail](../../../functions/mpi/sources/mpi_data_communicator/MPIDataCommunicator/GatherDetail.md)
- [GatherDetail](../../../functions/mpi/sources/mpi_data_communicator/MPIDataCommunicator/GatherDetail-2.md)
- [GathervDetail](../../../functions/mpi/sources/mpi_data_communicator/MPIDataCommunicator/GathervDetail.md)
- [GathervDetail](../../../functions/mpi/sources/mpi_data_communicator/MPIDataCommunicator/GathervDetail-2.md)
- [AllGatherDetail](../../../functions/mpi/sources/mpi_data_communicator/MPIDataCommunicator/AllGatherDetail.md)
- [AllGatherDetail](../../../functions/mpi/sources/mpi_data_communicator/MPIDataCommunicator/AllGatherDetail-2.md)
- [AllGathervDetail](../../../functions/mpi/sources/mpi_data_communicator/MPIDataCommunicator/AllGathervDetail.md)
- [AllGathervDetail](../../../functions/mpi/sources/mpi_data_communicator/MPIDataCommunicator/AllGathervDetail-2.md)
- [BroadcastErrorIfTrue](../../../functions/mpi/sources/mpi_data_communicator/MPIDataCommunicator/BroadcastErrorIfTrue.md)
- [BroadcastErrorIfFalse](../../../functions/mpi/sources/mpi_data_communicator/MPIDataCommunicator/BroadcastErrorIfFalse.md)
- [ErrorIfTrueOnAnyRank](../../../functions/mpi/sources/mpi_data_communicator/MPIDataCommunicator/ErrorIfTrueOnAnyRank.md)
- [ErrorIfFalseOnAnyRank](../../../functions/mpi/sources/mpi_data_communicator/MPIDataCommunicator/ErrorIfFalseOnAnyRank.md)
- [IsEqualOnAllRanks](../../../functions/mpi/sources/mpi_data_communicator/MPIDataCommunicator/IsEqualOnAllRanks.md)
- [IsValidRank](../../../functions/mpi/sources/mpi_data_communicator/MPIDataCommunicator/IsValidRank.md)
- [ValidateScattervInput](../../../functions/mpi/sources/mpi_data_communicator/MPIDataCommunicator/ValidateScattervInput.md)
- [ValidateGathervInput](../../../functions/mpi/sources/mpi_data_communicator/MPIDataCommunicator/ValidateGathervInput.md)
- [ValidateAllGathervInput](../../../functions/mpi/sources/mpi_data_communicator/MPIDataCommunicator/ValidateAllGathervInput.md)
- [PrepareScattervBuffers](../../../functions/mpi/sources/mpi_data_communicator/MPIDataCommunicator/PrepareScattervBuffers.md)
- [PrepareGathervBuffers](../../../functions/mpi/sources/mpi_data_communicator/MPIDataCommunicator/PrepareGathervBuffers.md)
- [PrepareAllGathervBuffers](../../../functions/mpi/sources/mpi_data_communicator/MPIDataCommunicator/PrepareAllGathervBuffers.md)
- [PrepareGathervReturn](../../../functions/mpi/sources/mpi_data_communicator/MPIDataCommunicator/PrepareGathervReturn.md)
- [PrepareAllGathervReturn](../../../functions/mpi/sources/mpi_data_communicator/MPIDataCommunicator/PrepareAllGathervReturn.md)
- [MPIDatatype](../../../functions/mpi/sources/mpi_data_communicator/MPIDataCommunicator/MPIDatatype.md)
- [MPIBuffer](../../../functions/mpi/sources/mpi_data_communicator/MPIDataCommunicator/MPIBuffer.md)
- [MPIBuffer](../../../functions/mpi/sources/mpi_data_communicator/MPIDataCommunicator/MPIBuffer-2.md)
- [MPIMessageSize](../../../functions/mpi/sources/mpi_data_communicator/MPIDataCommunicator/MPIMessageSize.md)

# Imports

- `<algorithm>`
- `"includes/parallel_environment.h"`
- `"mpi/includes/mpi_data_communicator.h"`
- `"mpi/includes/mpi_manager.h"`
- `"mpi/includes/mpi_message.h"`
- `"mpi/utilities/data_communicator_factory.h"`

# Member of

- [KratosMultiphysics](../../../packages/KratosMultiphysics.md)