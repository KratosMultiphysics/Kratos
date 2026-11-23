---
type: C++ Method
title: ReadModelPart
resource: sources/model_part_io.cpp#L831-L910
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/sources/model_part_io/ModelPartIO/ResetInput
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/containers/flags/IsNot
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/model_part_io/ModelPartIO/ReadModelPartDataBlock
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/model_part_io/ModelPartIO/SkipBlock
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/model_part_io/ModelPartIO/ReadPropertiesBlock
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/model_part_io/ModelPartIO/ReadNodalDataBlock
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/model_part_io/ModelPartIO/ReadElementalDataBlock
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/model_part_io/ModelPartIO/ReadConditionalDataBlock
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/model_part_io/ModelPartIO/ReadConstraintDataBlock
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/model_part_io/ModelPartIO/ReadCommunicatorDataBlock
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/model_part_io/ModelPartIO/ReadMeshBlock
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/model_part_io/ModelPartIO/ReadSubModelPartBlock
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/timer/Timer/Stop
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void ModelPartIO::ReadModelPart(ModelPart & rThisModelPart)`

# Calls

- [ResetInput](../../../../functions/sources/model_part_io/ModelPartIO/ResetInput.md)
- [IsNot](../../../../functions/containers/flags/IsNot.md)
- [ReadModelPartDataBlock](../../../../functions/sources/model_part_io/ModelPartIO/ReadModelPartDataBlock.md)
- [SkipBlock](../../../../functions/sources/model_part_io/ModelPartIO/SkipBlock.md)
- [ReadPropertiesBlock](../../../../functions/sources/model_part_io/ModelPartIO/ReadPropertiesBlock.md)
- [ReadNodalDataBlock](../../../../functions/sources/model_part_io/ModelPartIO/ReadNodalDataBlock.md)
- [ReadElementalDataBlock](../../../../functions/sources/model_part_io/ModelPartIO/ReadElementalDataBlock.md)
- [ReadConditionalDataBlock](../../../../functions/sources/model_part_io/ModelPartIO/ReadConditionalDataBlock.md)
- [ReadConstraintDataBlock](../../../../functions/sources/model_part_io/ModelPartIO/ReadConstraintDataBlock.md)
- [ReadCommunicatorDataBlock](../../../../functions/sources/model_part_io/ModelPartIO/ReadCommunicatorDataBlock.md)
- [ReadMeshBlock](../../../../functions/sources/model_part_io/ModelPartIO/ReadMeshBlock.md)
- [ReadSubModelPartBlock](../../../../functions/sources/model_part_io/ModelPartIO/ReadSubModelPartBlock.md)
- [Stop](../../../../functions/utilities/timer/Timer/Stop.md)