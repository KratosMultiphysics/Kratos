---
type: C++ Method
title: DivideMeshDataBlock
resource: sources/model_part_io.cpp#L5697-L5709
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/sources/model_part_io/ModelPartIO/WriteInAllFiles
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/sources/model_part_io/ModelPartIO/DivideMeshBlock
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void ModelPartIO::DivideMeshDataBlock(OutputFilesContainerType& OutputFiles)`

# Calls

- [WriteInAllFiles](../../../../functions/sources/model_part_io/ModelPartIO/WriteInAllFiles.md)

# Called by

- [DivideMeshBlock](../../../../functions/sources/model_part_io/ModelPartIO/DivideMeshBlock.md)