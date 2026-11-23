---
type: C++ Method
title: RemoveElementAndBelongings
resource: utilities/auxiliar_model_part_utilities.cpp#L193-L224
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
  - target: functions/sources/model_part/ModelPart/RemoveConditions
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/model_part/ModelPart/RemoveNodes
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void AuxiliarModelPartUtilities::RemoveElementAndBelongings( IndexType ElementId, Flags IdentifierFlag, IndexType ThisIndex)`

# Calls

- [SetFlag](../../../../functions/utilities/variable_utils/SetFlag.md)
- [GetGeometry](../../../../functions/input_output/cad_json_input/CadJsonInput/GetGeometry.md)
- [IsNot](../../../../functions/containers/flags/IsNot.md)
- [RemoveConditions](../../../../functions/sources/model_part/ModelPart/RemoveConditions.md)
- [RemoveNodes](../../../../functions/sources/model_part/ModelPart/RemoveNodes.md)