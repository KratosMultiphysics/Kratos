---
type: C++ Module
title: data_communicator
resource: includes/data_communicator.h#L1-L1060
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
  - target: external/type-traits
    resolved_by: tree-sitter
    confidence: exact
  - target: external/containers-array-1d-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/containers-flags-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/includes-define-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/includes-mpi-serializer-h
    resolved_by: tree-sitter
    confidence: exact
  member_of:
  - target: packages/KratosMultiphysics
    resolved_by: tree-sitter
    confidence: exact
---

# Contains

- [TypeFromBool](../../classes/includes/data_communicator/TypeFromBool.md)
- [CheckSerializationForSimpleType](../../functions/includes/data_communicator/CheckSerializationForSimpleType.md)
- [CheckSerializationForSimpleType](../../functions/includes/data_communicator/CheckSerializationForSimpleType-2.md)
- [DataCommunicator](../../functions/includes/data_communicator/DataCommunicator.md)
- [Create](../../functions/includes/data_communicator/Create.md)
- [Barrier](../../functions/includes/data_communicator/Barrier.md)
- [AndReduce](../../functions/includes/data_communicator/AndReduce.md)
- [AndReduce](../../functions/includes/data_communicator/AndReduce-2.md)
- [OrReduce](../../functions/includes/data_communicator/OrReduce.md)
- [OrReduce](../../functions/includes/data_communicator/OrReduce-2.md)
- [AndReduceAll](../../functions/includes/data_communicator/AndReduceAll.md)
- [AndReduceAll](../../functions/includes/data_communicator/AndReduceAll-2.md)
- [OrReduceAll](../../functions/includes/data_communicator/OrReduceAll.md)
- [OrReduceAll](../../functions/includes/data_communicator/OrReduceAll-2.md)
- [Broadcast](../../functions/includes/data_communicator/Broadcast.md)
- [SendRecv](../../functions/includes/data_communicator/SendRecv.md)
- [SendRecv](../../functions/includes/data_communicator/SendRecv-2.md)
- [SendRecv](../../functions/includes/data_communicator/SendRecv-3.md)
- [SendRecv](../../functions/includes/data_communicator/SendRecv-4.md)
- [Send](../../functions/includes/data_communicator/Send.md)
- [Recv](../../functions/includes/data_communicator/Recv.md)
- [Rank](../../functions/includes/data_communicator/Rank.md)
- [Size](../../functions/includes/data_communicator/Size.md)
- [IsDistributed](../../functions/includes/data_communicator/IsDistributed.md)
- [IsDefinedOnThisRank](../../functions/includes/data_communicator/IsDefinedOnThisRank.md)
- [IsNullOnThisRank](../../functions/includes/data_communicator/IsNullOnThisRank.md)
- [BroadcastErrorIfTrue](../../functions/includes/data_communicator/BroadcastErrorIfTrue.md)
- [BroadcastErrorIfFalse](../../functions/includes/data_communicator/BroadcastErrorIfFalse.md)
- [ErrorIfTrueOnAnyRank](../../functions/includes/data_communicator/ErrorIfTrueOnAnyRank.md)
- [ErrorIfFalseOnAnyRank](../../functions/includes/data_communicator/ErrorIfFalseOnAnyRank.md)
- [Info](../../functions/includes/data_communicator/Info.md)
- [PrintInfo](../../functions/includes/data_communicator/PrintInfo.md)
- [PrintData](../../functions/includes/data_communicator/PrintData.md)
- [BroadcastImpl](../../functions/includes/data_communicator/BroadcastImpl.md)
- [SendRecvImpl](../../functions/includes/data_communicator/SendRecvImpl.md)
- [SendRecvImpl](../../functions/includes/data_communicator/SendRecvImpl-2.md)
- [SendRecvImpl](../../functions/includes/data_communicator/SendRecvImpl-3.md)
- [SendImpl](../../functions/includes/data_communicator/SendImpl.md)
- [SendImpl](../../functions/includes/data_communicator/SendImpl-2.md)
- [RecvImpl](../../functions/includes/data_communicator/RecvImpl.md)
- [RecvImpl](../../functions/includes/data_communicator/RecvImpl-2.md)

# Imports

- `<string>`
- `<iostream>`
- `<type_traits>`
- `"containers/array_1d.h"`
- `"containers/flags.h"`
- `"includes/define.h"`
- `"includes/mpi_serializer.h"`

# Member of

- [KratosMultiphysics](../../packages/KratosMultiphysics.md)