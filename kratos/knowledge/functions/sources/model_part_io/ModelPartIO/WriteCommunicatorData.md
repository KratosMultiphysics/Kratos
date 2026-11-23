---
type: C++ Method
title: WriteCommunicatorData
resource: sources/model_part_io.cpp#L6201-L6333
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/sources/model_part_io/ModelPartIO/WriteInAllFiles
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/sources/model_part_io/ModelPartIO/DivideInputToPartitionsImpl
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void ModelPartIO::WriteCommunicatorData( OutputFilesContainerType& OutputFiles, SizeType NumberOfPartitions, GraphType const& DomainsColoredGraph, PartitionIndicesType const& NodesPartitions, PartitionIndicesType const& ElementsPartitions, PartitionIndicesType const& ConditionsPartitions, PartitionIndicesType const& ConstraintsPartitions, PartitionIndicesType const& GeometriesPartitions, PartitionIndicesContainerType const& NodesAllPartitions, PartitionIndicesContainerType const& ElementsAllPartitions, PartitionIndicesContainerType const& ConditionsAllPartitions, PartitionIndicesContainerType const& rConstraintsAllPartitions, PartitionIndicesContainerType const& GeometriesAllPartitions )`

# Calls

- [WriteInAllFiles](../../../../functions/sources/model_part_io/ModelPartIO/WriteInAllFiles.md)

# Called by

- [DivideInputToPartitionsImpl](../../../../functions/sources/model_part_io/ModelPartIO/DivideInputToPartitionsImpl.md)