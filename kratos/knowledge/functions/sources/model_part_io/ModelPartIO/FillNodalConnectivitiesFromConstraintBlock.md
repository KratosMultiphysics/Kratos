---
type: C++ Method
title: FillNodalConnectivitiesFromConstraintBlock
resource: sources/model_part_io.cpp#L3839-L3941
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
  called_by:
  - target: functions/sources/model_part_io/ModelPartIO/ReadNodalGraph
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void ModelPartIO::FillNodalConnectivitiesFromConstraintBlock(ConnectivitiesContainerType& rNodalConnectivities)`

# Calls

- [SplitStringIntoAVector](../../../../functions/utilities/string_utilities/SplitStringIntoAVector.md)
- [CountValuesUntilCharacter](../../../../functions/utilities/string_utilities/CountValuesUntilCharacter.md)
- [CheckEndBlock](../../../../functions/sources/model_part_io/ModelPartIO/CheckEndBlock.md)

# Called by

- [ReadNodalGraph](../../../../functions/sources/model_part_io/ModelPartIO/ReadNodalGraph.md)