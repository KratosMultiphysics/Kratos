---
type: C++ Function
title: norm_frobenius
resource: includes/amatrix_interface.h#L687-L696
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/includes/mortar_classes/DualLagrangeMultiplierOperators/CalculateAe
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/includes/mortar_classes/DualLagrangeMultiplierOperatorsWithDerivatives/CalculateDeltaAe
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/simple_mortar_mapper_process/PerformMortarOperations
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/spaces/kratos_space/KratosSpace/TwoNorm
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/geometries/test_geometry/VerifyStrainExactness
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/geometry_tester/GeometryTesterUtility/VerifyStrainExactness
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/geometry_tester/GeometryTesterUtility/VerifyShapeFunctionsSecondDerivativesValues
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/math_utils/CheckConditionNumber
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/math_utils/CalculateExponentialOfMatrix
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`typename TExpressionType::data_type norm_frobenius( AMatrix::MatrixExpression<TExpressionType, TCategory> const& TheExpression)`

# Called by

- [CalculateAe](../../../functions/includes/mortar_classes/DualLagrangeMultiplierOperators/CalculateAe.md)
- [CalculateDeltaAe](../../../functions/includes/mortar_classes/DualLagrangeMultiplierOperatorsWithDerivatives/CalculateDeltaAe.md)
- [PerformMortarOperations](../../../functions/processes/simple_mortar_mapper_process/PerformMortarOperations.md)
- [TwoNorm](../../../functions/spaces/kratos_space/KratosSpace/TwoNorm.md)
- [VerifyStrainExactness](../../../functions/tests/cpp_tests/geometries/test_geometry/VerifyStrainExactness.md)
- [VerifyStrainExactness](../../../functions/utilities/geometry_tester/GeometryTesterUtility/VerifyStrainExactness.md)
- [VerifyShapeFunctionsSecondDerivativesValues](../../../functions/utilities/geometry_tester/GeometryTesterUtility/VerifyShapeFunctionsSecondDerivativesValues.md)
- [CheckConditionNumber](../../../functions/utilities/math_utils/CheckConditionNumber.md)
- [CalculateExponentialOfMatrix](../../../functions/utilities/math_utils/CalculateExponentialOfMatrix.md)