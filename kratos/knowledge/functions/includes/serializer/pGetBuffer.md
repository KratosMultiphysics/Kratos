---
type: C++ Function
title: pGetBuffer
resource: includes/serializer.h#L960-L963
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/includes/stream_serializer/StreamSerializer/StreamSerializer
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/includes/stream_serializer/StreamSerializer/GetStringRepresentation
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/mpi/includes/mpi_communicator/AsyncSendAndReceiveObjects
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/mpi/includes/mpi_communicator/UpdateValues
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/mpi/utilities/gather_modelpart_utility/GatherModelPartUtility/GatherEntityFromOtherPartitions
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python/add_serializer_to_python/SerializerPrint
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/file_serializer/FileSerializer/FileSerializer
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`BufferType* pGetBuffer()`

# Called by

- [StreamSerializer](../../../functions/includes/stream_serializer/StreamSerializer/StreamSerializer.md)
- [GetStringRepresentation](../../../functions/includes/stream_serializer/StreamSerializer/GetStringRepresentation.md)
- [AsyncSendAndReceiveObjects](../../../functions/mpi/includes/mpi_communicator/AsyncSendAndReceiveObjects.md)
- [UpdateValues](../../../functions/mpi/includes/mpi_communicator/UpdateValues.md)
- [GatherEntityFromOtherPartitions](../../../functions/mpi/utilities/gather_modelpart_utility/GatherModelPartUtility/GatherEntityFromOtherPartitions.md)
- [SerializerPrint](../../../functions/python/add_serializer_to_python/SerializerPrint.md)
- [FileSerializer](../../../functions/sources/file_serializer/FileSerializer/FileSerializer.md)