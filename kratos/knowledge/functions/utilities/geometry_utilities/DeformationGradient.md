---
type: C++ Function
title: DeformationGradient
resource: utilities/geometry_utilities.h#L690-L700
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/includes/amatrix_interface/prod
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/tests/cpp_tests/utilities/test_geometry_utils/KRATOS_TEST_CASE_IN_SUITE
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`static void DeformationGradient( TMatrix1 const& rJ, TMatrix2 const& rInvJ0, TMatrix3& rF )`

# Calls

- [prod](../../../functions/includes/amatrix_interface/prod.md)

# Called by

- [KRATOS_TEST_CASE_IN_SUITE](../../../functions/tests/cpp_tests/utilities/test_geometry_utils/KRATOS_TEST_CASE_IN_SUITE.md)