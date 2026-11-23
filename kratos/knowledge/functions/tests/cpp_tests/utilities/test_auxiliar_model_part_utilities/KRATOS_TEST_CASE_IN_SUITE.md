---
type: C++ Function
title: KRATOS_TEST_CASE_IN_SUITE
resource: tests/cpp_tests/utilities/test_auxiliar_model_part_utilities.cpp#L233-L281
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/utilities/auxiliar_model_part_utilities/AuxiliarModelPartUtilities
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/auxiliar_model_part_utilities/AuxiliarModelPartUtilities/RemoveOrphanNodesFromSubModelParts
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/model_part/ModelPart/HasSubModelPart
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/model_part/ModelPart/CreateNewProperties
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/utilities/test_auxiliar_model_part_utilities/ComputationGetData
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/auxiliar_model_part_utilities/GetScalarData
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/auxiliar_model_part_utilities/GetVectorData
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/utilities/test_auxiliar_model_part_utilities/PreComputeSetData
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/auxiliar_model_part_utilities/SetScalarData
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/utilities/test_auxiliar_model_part_utilities/PostComputeSetData
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/auxiliar_model_part_utilities/SetVectorData
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`KRATOS_TEST_CASE_IN_SUITE(AuxiliarModelPartUtilities_RemoveNodesFromSubModePartsWithoutCorrespondingEntities, KratosCoreFastSuite)`

# Calls

- [AuxiliarModelPartUtilities](../../../../../functions/utilities/auxiliar_model_part_utilities/AuxiliarModelPartUtilities.md)
- [RemoveOrphanNodesFromSubModelParts](../../../../../functions/utilities/auxiliar_model_part_utilities/AuxiliarModelPartUtilities/RemoveOrphanNodesFromSubModelParts.md)
- [HasSubModelPart](../../../../../functions/sources/model_part/ModelPart/HasSubModelPart.md)
- [CreateNewProperties](../../../../../functions/sources/model_part/ModelPart/CreateNewProperties.md)
- [ComputationGetData](../../../../../functions/tests/cpp_tests/utilities/test_auxiliar_model_part_utilities/ComputationGetData.md)
- [GetScalarData](../../../../../functions/utilities/auxiliar_model_part_utilities/GetScalarData.md)
- [GetVectorData](../../../../../functions/utilities/auxiliar_model_part_utilities/GetVectorData.md)
- [PreComputeSetData](../../../../../functions/tests/cpp_tests/utilities/test_auxiliar_model_part_utilities/PreComputeSetData.md)
- [SetScalarData](../../../../../functions/utilities/auxiliar_model_part_utilities/SetScalarData.md)
- [PostComputeSetData](../../../../../functions/tests/cpp_tests/utilities/test_auxiliar_model_part_utilities/PostComputeSetData.md)
- [SetVectorData](../../../../../functions/utilities/auxiliar_model_part_utilities/SetVectorData.md)