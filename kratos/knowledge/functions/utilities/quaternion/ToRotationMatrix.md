---
type: C++ Function
title: ToRotationMatrix
resource: utilities/quaternion.h#L213-L229
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/utilities/quaternion/W-2
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/tests/cpp_tests/utilities/test_quaternion/KRATOS_TEST_CASE_IN_SUITE
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/beam_math_utilities/ExponentialTransform
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/beam_math_utilities/BeamMathUtils/MapLocalToGlobal2D
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/beam_math_utilities/BeamMathUtils/MapLocalToGlobal3D
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`inline void ToRotationMatrix(TMatrix3x3& R)const`

# Calls

- [W](../../../functions/utilities/quaternion/W-2.md)

# Called by

- [KRATOS_TEST_CASE_IN_SUITE](../../../functions/tests/cpp_tests/utilities/test_quaternion/KRATOS_TEST_CASE_IN_SUITE.md)
- [ExponentialTransform](../../../functions/utilities/beam_math_utilities/ExponentialTransform.md)
- [MapLocalToGlobal2D](../../../functions/utilities/beam_math_utilities/BeamMathUtils/MapLocalToGlobal2D.md)
- [MapLocalToGlobal3D](../../../functions/utilities/beam_math_utilities/BeamMathUtils/MapLocalToGlobal3D.md)