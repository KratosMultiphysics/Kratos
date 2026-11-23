---
type: C++ Method
title: DivideInputToPartitions
resource: sources/model_part_io.cpp#L1415-L1450
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/sources/kratos_filesystem/parent_path
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kratos_filesystem/FilesystemExtensions/MPISafeCreateDirectories
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/model_part_io/ModelPartIO/DivideInputToPartitionsImpl
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void ModelPartIO::DivideInputToPartitions( SizeType NumberOfPartitions, const PartitioningInfo& rPartitioningInfo)`

# Calls

- [parent_path](../../../../functions/sources/kratos_filesystem/parent_path.md)
- [MPISafeCreateDirectories](../../../../functions/sources/kratos_filesystem/FilesystemExtensions/MPISafeCreateDirectories.md)
- [DivideInputToPartitionsImpl](../../../../functions/sources/model_part_io/ModelPartIO/DivideInputToPartitionsImpl.md)