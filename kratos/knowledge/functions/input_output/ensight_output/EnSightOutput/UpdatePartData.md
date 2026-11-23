---
type: C++ Method
title: UpdatePartData
resource: input_output/ensight_output.cpp#L2503-L2564
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/sources/kratos_parameters/Parameters/GetBool
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/input_output/ensight_output/EnSightOutput/CollectPartData
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/input_output/ensight_output/EnSightOutput/EnSightOutput
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/input_output/ensight_output/EnSightOutput/PrintOutput
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void EnSightOutput::UpdatePartData()`

# Calls

- [GetBool](../../../../functions/sources/kratos_parameters/Parameters/GetBool.md)
- [CollectPartData](../../../../functions/input_output/ensight_output/EnSightOutput/CollectPartData.md)

# Called by

- [EnSightOutput](../../../../functions/input_output/ensight_output/EnSightOutput/EnSightOutput.md)
- [PrintOutput](../../../../functions/input_output/ensight_output/EnSightOutput/PrintOutput.md)