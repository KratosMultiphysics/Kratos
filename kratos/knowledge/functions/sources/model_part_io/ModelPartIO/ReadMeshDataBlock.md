---
type: C++ Method
title: ReadMeshDataBlock
resource: sources/model_part_io.cpp#L4181-L4256
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/sources/model_part_io/ModelPartIO/CheckEndBlock
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/sources/model_part_io/ModelPartIO/ReadMeshBlock
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void ModelPartIO::ReadMeshDataBlock(MeshType& rMesh)`

# Calls

- [CheckEndBlock](../../../../functions/sources/model_part_io/ModelPartIO/CheckEndBlock.md)

# Called by

- [ReadMeshBlock](../../../../functions/sources/model_part_io/ModelPartIO/ReadMeshBlock.md)