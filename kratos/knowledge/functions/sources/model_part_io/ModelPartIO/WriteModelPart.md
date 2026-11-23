---
type: C++ Method
title: WriteModelPart
resource: sources/model_part_io.cpp#L912-L943
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/utilities/timer/Timer/Start
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/containers/flags/IsNot
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/model_part_io/ModelPartIO/WriteNodalDataBlock
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/model_part_io/ModelPartIO/WriteSubModelPartBlock
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/timer/Timer/Stop
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void ModelPartIO::WriteModelPart(ModelPart& rThisModelPart)`

# Calls

- [Start](../../../../functions/utilities/timer/Timer/Start.md)
- [IsNot](../../../../functions/containers/flags/IsNot.md)
- [WriteNodalDataBlock](../../../../functions/sources/model_part_io/ModelPartIO/WriteNodalDataBlock.md)
- [WriteSubModelPartBlock](../../../../functions/sources/model_part_io/ModelPartIO/WriteSubModelPartBlock.md)
- [Stop](../../../../functions/utilities/timer/Timer/Stop.md)