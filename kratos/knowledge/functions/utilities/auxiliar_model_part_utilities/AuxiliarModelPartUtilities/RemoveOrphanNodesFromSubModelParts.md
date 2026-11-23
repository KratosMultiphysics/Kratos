---
type: C++ Method
title: RemoveOrphanNodesFromSubModelParts
resource: utilities/auxiliar_model_part_utilities.cpp#L453-L485
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
  - target: functions/sources/model_part/ModelPart/RemoveNodes
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/tests/cpp_tests/utilities/test_auxiliar_model_part_utilities/KRATOS_TEST_CASE_IN_SUITE
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void AuxiliarModelPartUtilities::RemoveOrphanNodesFromSubModelParts()`

# Calls

- [SetFlag](../../../../functions/utilities/variable_utils/SetFlag.md)
- [GetGeometry](../../../../functions/input_output/cad_json_input/CadJsonInput/GetGeometry.md)
- [RemoveNodes](../../../../functions/sources/model_part/ModelPart/RemoveNodes.md)

# Called by

- [KRATOS_TEST_CASE_IN_SUITE](../../../../functions/tests/cpp_tests/utilities/test_auxiliar_model_part_utilities/KRATOS_TEST_CASE_IN_SUITE.md)