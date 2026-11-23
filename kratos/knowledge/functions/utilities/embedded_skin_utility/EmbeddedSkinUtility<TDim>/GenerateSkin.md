---
type: C++ Method
title: GenerateSkin
resource: utilities/embedded_skin_utility.cpp#L30-L64
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/utilities/embedded_skin_utility/EmbeddedSkinUtility<TDim>/SetSkinEntitiesProperties
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/input_output/cad_json_input/CadJsonInput/GetGeometry
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/embedded_skin_utility/EmbeddedSkinUtility<TDim>/SetDistancesVector
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/embedded_skin_utility/EmbeddedSkinUtility<TDim>/ElementIsSplit
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/embedded_skin_utility/EmbeddedSkinUtility<TDim>/ComputeElementSkin
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/embedded_skin_utility/EmbeddedSkinUtility<TDim>/RenumberAndAddSkinEntities
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/tests/cpp_tests/utilities/test_embedded_skin_utility/KRATOS_TEST_CASE_IN_SUITE
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_embedded_skin_mapping/TestEmbeddedSkinMapping/test_embedded_skin_mapping
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void EmbeddedSkinUtility<TDim>::GenerateSkin()`

# Calls

- [SetSkinEntitiesProperties](../../../../functions/utilities/embedded_skin_utility/EmbeddedSkinUtility<TDim>/SetSkinEntitiesProperties.md)
- [GetGeometry](../../../../functions/input_output/cad_json_input/CadJsonInput/GetGeometry.md)
- [SetDistancesVector](../../../../functions/utilities/embedded_skin_utility/EmbeddedSkinUtility<TDim>/SetDistancesVector.md)
- [ElementIsSplit](../../../../functions/utilities/embedded_skin_utility/EmbeddedSkinUtility<TDim>/ElementIsSplit.md)
- [ComputeElementSkin](../../../../functions/utilities/embedded_skin_utility/EmbeddedSkinUtility<TDim>/ComputeElementSkin.md)
- [RenumberAndAddSkinEntities](../../../../functions/utilities/embedded_skin_utility/EmbeddedSkinUtility<TDim>/RenumberAndAddSkinEntities.md)

# Called by

- [KRATOS_TEST_CASE_IN_SUITE](../../../../functions/tests/cpp_tests/utilities/test_embedded_skin_utility/KRATOS_TEST_CASE_IN_SUITE.md)
- [test_embedded_skin_mapping](../../../../functions/tests/test_embedded_skin_mapping/TestEmbeddedSkinMapping/test_embedded_skin_mapping.md)