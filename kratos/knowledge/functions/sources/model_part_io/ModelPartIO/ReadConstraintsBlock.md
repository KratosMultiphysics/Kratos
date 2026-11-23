---
type: C++ Method
title: ReadConstraintsBlock
resource: sources/model_part_io.cpp#L2410-L2419
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/utilities/string_utilities/SplitStringIntoAVector
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/string_utilities/CountValuesUntilCharacter
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/model_part_io/ModelPartIO/CheckEndBlock
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/model_part_io/ModelPartIO/FindKey
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/includes/node/HasDofFor
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void ModelPartIO::ReadConstraintsBlock(ModelPart& rModelPart)`

# Calls

- [SplitStringIntoAVector](../../../../functions/utilities/string_utilities/SplitStringIntoAVector.md)
- [CountValuesUntilCharacter](../../../../functions/utilities/string_utilities/CountValuesUntilCharacter.md)
- [CheckEndBlock](../../../../functions/sources/model_part_io/ModelPartIO/CheckEndBlock.md)
- [FindKey](../../../../functions/sources/model_part_io/ModelPartIO/FindKey.md)
- [HasDofFor](../../../../functions/includes/node/HasDofFor.md)