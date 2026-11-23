---
type: C++ Method
title: RemoveConditionsAndBelongings
resource: utilities/auxiliar_model_part_utilities.cpp#L403-L438
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/utilities/variable_utils/SetFlag
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/containers/flags/IsNot
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/input_output/cad_json_input/CadJsonInput/GetGeometry
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/model_part/ModelPart/RemoveConditions
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/model_part/ModelPart/RemoveElements
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/model_part/ModelPart/RemoveNodes
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/tests/cpp_tests/sources/test_model_part/KRATOS_TEST_CASE_IN_SUITE
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/auxiliar_model_part_utilities/AuxiliarModelPartUtilities/RemoveConditionsAndBelongingsFromAllLevels
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void AuxiliarModelPartUtilities::RemoveConditionsAndBelongings(Flags IdentifierFlag)`

# Calls

- [SetFlag](../../../../functions/utilities/variable_utils/SetFlag.md)
- [IsNot](../../../../functions/containers/flags/IsNot.md)
- [GetGeometry](../../../../functions/input_output/cad_json_input/CadJsonInput/GetGeometry.md)
- [RemoveConditions](../../../../functions/sources/model_part/ModelPart/RemoveConditions.md)
- [RemoveElements](../../../../functions/sources/model_part/ModelPart/RemoveElements.md)
- [RemoveNodes](../../../../functions/sources/model_part/ModelPart/RemoveNodes.md)

# Called by

- [KRATOS_TEST_CASE_IN_SUITE](../../../../functions/tests/cpp_tests/sources/test_model_part/KRATOS_TEST_CASE_IN_SUITE.md)
- [RemoveConditionsAndBelongingsFromAllLevels](../../../../functions/utilities/auxiliar_model_part_utilities/AuxiliarModelPartUtilities/RemoveConditionsAndBelongingsFromAllLevels.md)