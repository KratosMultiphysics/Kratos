---
type: C++ Method
title: RemoveConditionAndBelongings
resource: utilities/auxiliar_model_part_utilities.cpp#L324-L353
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/utilities/variable_utils/SetFlag
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/input_output/cad_json_input/CadJsonInput/GetGeometry
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/containers/flags/IsNot
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/model_part/ModelPart/RemoveElements
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/model_part/ModelPart/RemoveNodes
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void AuxiliarModelPartUtilities::RemoveConditionAndBelongings(IndexType ConditionId, Flags IdentifierFlag, IndexType ThisIndex)`

# Calls

- [SetFlag](../../../../functions/utilities/variable_utils/SetFlag.md)
- [GetGeometry](../../../../functions/input_output/cad_json_input/CadJsonInput/GetGeometry.md)
- [IsNot](../../../../functions/containers/flags/IsNot.md)
- [RemoveElements](../../../../functions/sources/model_part/ModelPart/RemoveElements.md)
- [RemoveNodes](../../../../functions/sources/model_part/ModelPart/RemoveNodes.md)