---
type: C++ Function
title: CofactorMatrix
resource: utilities/math_utils.h#L191-L202
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/utilities/math_utils/Cofactor
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/tests/cpp_tests/utilities/test_math_utils/KRATOS_TEST_CASE_IN_SUITE
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/geometrical_sensitivity_utility/GeometricalSensitivityUtility/Initialize
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/geometrical_sensitivity_utility/GeometricalSensitivityUtility/CalculateCofactorOfJacobianSensitivity
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/line_sensitivity_utility/LineSensitivityUtility/Initialize
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`static MatrixType CofactorMatrix(const TMatrixType& rMat)`

# Calls

- [Cofactor](../../../functions/utilities/math_utils/Cofactor.md)

# Called by

- [KRATOS_TEST_CASE_IN_SUITE](../../../functions/tests/cpp_tests/utilities/test_math_utils/KRATOS_TEST_CASE_IN_SUITE.md)
- [Initialize](../../../functions/utilities/geometrical_sensitivity_utility/GeometricalSensitivityUtility/Initialize.md)
- [CalculateCofactorOfJacobianSensitivity](../../../functions/utilities/geometrical_sensitivity_utility/GeometricalSensitivityUtility/CalculateCofactorOfJacobianSensitivity.md)
- [Initialize](../../../functions/utilities/line_sensitivity_utility/LineSensitivityUtility/Initialize.md)