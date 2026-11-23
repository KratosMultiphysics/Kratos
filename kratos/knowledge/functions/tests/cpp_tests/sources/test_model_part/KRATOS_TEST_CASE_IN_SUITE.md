---
type: C++ Function
title: KRATOS_TEST_CASE_IN_SUITE
resource: tests/cpp_tests/sources/test_model_part.cpp#L58-L66
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/includes/model_part/NumberOfSubModelParts
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/model_part/ModelPart/GetSubModelPartNames
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/includes/model_part/FullName
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/sources/test_model_part/GenerateGenericModelPart
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/model_part/ModelPart/AddTable
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/includes/table/Table<double, double>/GetDerivative
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/model_part/ModelPart/RemoveElements
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/auxiliar_model_part_utilities/AuxiliarModelPartUtilities
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/auxiliar_model_part_utilities/AuxiliarModelPartUtilities/RemoveElementsAndBelongings
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/model_part/ModelPart/RemoveConditions
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/auxiliar_model_part_utilities/AuxiliarModelPartUtilities/RemoveConditionsAndBelongings
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/auxiliar_model_part_utilities/AuxiliarModelPartUtilities/EnsureModelPartOwnsProperties
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/model_part/ModelPart/HasSubModelPart
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/model_part/ModelPart/CreateNewProperties
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`KRATOS_TEST_CASE_IN_SUITE(ModelPartDataValueContainer, KratosCoreFastSuite)`

# Calls

- [NumberOfSubModelParts](../../../../../functions/includes/model_part/NumberOfSubModelParts.md)
- [GetSubModelPartNames](../../../../../functions/sources/model_part/ModelPart/GetSubModelPartNames.md)
- [FullName](../../../../../functions/includes/model_part/FullName.md)
- [GenerateGenericModelPart](../../../../../functions/tests/cpp_tests/sources/test_model_part/GenerateGenericModelPart.md)
- [AddTable](../../../../../functions/sources/model_part/ModelPart/AddTable.md)
- [GetDerivative](../../../../../functions/includes/table/Table<double, double>/GetDerivative.md)
- [RemoveElements](../../../../../functions/sources/model_part/ModelPart/RemoveElements.md)
- [AuxiliarModelPartUtilities](../../../../../functions/utilities/auxiliar_model_part_utilities/AuxiliarModelPartUtilities.md)
- [RemoveElementsAndBelongings](../../../../../functions/utilities/auxiliar_model_part_utilities/AuxiliarModelPartUtilities/RemoveElementsAndBelongings.md)
- [RemoveConditions](../../../../../functions/sources/model_part/ModelPart/RemoveConditions.md)
- [RemoveConditionsAndBelongings](../../../../../functions/utilities/auxiliar_model_part_utilities/AuxiliarModelPartUtilities/RemoveConditionsAndBelongings.md)
- [EnsureModelPartOwnsProperties](../../../../../functions/utilities/auxiliar_model_part_utilities/AuxiliarModelPartUtilities/EnsureModelPartOwnsProperties.md)
- [HasSubModelPart](../../../../../functions/sources/model_part/ModelPart/HasSubModelPart.md)
- [CreateNewProperties](../../../../../functions/sources/model_part/ModelPart/CreateNewProperties.md)