---
type: C++ Method
title: ReadConstraintsConnectivitiesBlock
resource: sources/model_part_io.cpp#L3532-L3642
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
  - target: functions/sources/model_part_io/ModelPartIO/ReorderedConstraintId
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/sources/model_part_io/ModelPartIO/ReadConstraintsConnectivities
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`ModelPartIO::SizeType ModelPartIO::ReadConstraintsConnectivitiesBlock(ConnectivitiesContainerType& rThisConnectivities)`

# Calls

- [SplitStringIntoAVector](../../../../functions/utilities/string_utilities/SplitStringIntoAVector.md)
- [CountValuesUntilCharacter](../../../../functions/utilities/string_utilities/CountValuesUntilCharacter.md)
- [CheckEndBlock](../../../../functions/sources/model_part_io/ModelPartIO/CheckEndBlock.md)
- [ReorderedConstraintId](../../../../functions/sources/model_part_io/ModelPartIO/ReorderedConstraintId.md)

# Called by

- [ReadConstraintsConnectivities](../../../../functions/sources/model_part_io/ModelPartIO/ReadConstraintsConnectivities.md)