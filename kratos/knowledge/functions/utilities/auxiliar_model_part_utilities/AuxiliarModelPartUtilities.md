---
type: C++ Function
title: AuxiliarModelPartUtilities
resource: utilities/auxiliar_model_part_utilities.h#L70-L73
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/processes/check_same_modelpart_using_skin_distance_process/CheckSameModelPartUsingSkinDistanceProcess<TDim>/Execute
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/sources/test_model_part/KRATOS_TEST_CASE_IN_SUITE
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/utilities/test_auxiliar_model_part_utilities/KRATOS_TEST_CASE_IN_SUITE
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_find_conditions_neighbours_process/TestFindConditionsNeighboursProcess/test_sphere
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/auxiliar_model_part_utilities/AuxiliarModelPartUtilities/RecursiveEnsureModelPartOwnsProperties
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/auxiliar_model_part_utilities/AuxiliarModelPartUtilities/RemoveElementAndBelongingsFromAllLevels
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/auxiliar_model_part_utilities/AuxiliarModelPartUtilities/RemoveElementsAndBelongingsFromAllLevels
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/auxiliar_model_part_utilities/AuxiliarModelPartUtilities/RemoveConditionAndBelongingsFromAllLevels
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/auxiliar_model_part_utilities/AuxiliarModelPartUtilities/RemoveConditionsAndBelongingsFromAllLevels
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`AuxiliarModelPartUtilities(ModelPart& rModelPart): mrModelPart(rModelPart)`

# Called by

- [Execute](../../../functions/processes/check_same_modelpart_using_skin_distance_process/CheckSameModelPartUsingSkinDistanceProcess<TDim>/Execute.md)
- [KRATOS_TEST_CASE_IN_SUITE](../../../functions/tests/cpp_tests/sources/test_model_part/KRATOS_TEST_CASE_IN_SUITE.md)
- [KRATOS_TEST_CASE_IN_SUITE](../../../functions/tests/cpp_tests/utilities/test_auxiliar_model_part_utilities/KRATOS_TEST_CASE_IN_SUITE.md)
- [test_sphere](../../../functions/tests/test_find_conditions_neighbours_process/TestFindConditionsNeighboursProcess/test_sphere.md)
- [RecursiveEnsureModelPartOwnsProperties](../../../functions/utilities/auxiliar_model_part_utilities/AuxiliarModelPartUtilities/RecursiveEnsureModelPartOwnsProperties.md)
- [RemoveElementAndBelongingsFromAllLevels](../../../functions/utilities/auxiliar_model_part_utilities/AuxiliarModelPartUtilities/RemoveElementAndBelongingsFromAllLevels.md)
- [RemoveElementsAndBelongingsFromAllLevels](../../../functions/utilities/auxiliar_model_part_utilities/AuxiliarModelPartUtilities/RemoveElementsAndBelongingsFromAllLevels.md)
- [RemoveConditionAndBelongingsFromAllLevels](../../../functions/utilities/auxiliar_model_part_utilities/AuxiliarModelPartUtilities/RemoveConditionAndBelongingsFromAllLevels.md)
- [RemoveConditionsAndBelongingsFromAllLevels](../../../functions/utilities/auxiliar_model_part_utilities/AuxiliarModelPartUtilities/RemoveConditionsAndBelongingsFromAllLevels.md)