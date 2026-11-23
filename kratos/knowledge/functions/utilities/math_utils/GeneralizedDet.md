---
type: C++ Function
title: GeneralizedDet
resource: utilities/math_utils.h#L650-L661
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/includes/amatrix_interface/prod
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/geometries/geometry/DeterminantOfJacobian
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/geometries/nurbs_surface_on_volume_geometry/NurbsSurfaceOnVolumeGeometry/DeterminantOfJacobian
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/geometries/quadrature_point_surface_in_volume_geometry/QuadraturePointSurfaceInVolumeGeometry/DeterminantOfJacobian
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/geometries/surface_in_nurbs_volume_geometry/SurfaceInNurbsVolumeGeometry/DeterminantOfJacobian
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/calculate_nodal_area_process/CalculateNodalAreaProcess<THistorical>/Execute
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/geometries/test_non_rectangular_jacobian/KRATOS_TEST_CASE_IN_SUITE
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/utilities/test_math_utils/KRATOS_TEST_CASE_IN_SUITE
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`static inline double GeneralizedDet(const TMatrixType& rA)`

# Calls

- [prod](../../../functions/includes/amatrix_interface/prod.md)

# Called by

- [DeterminantOfJacobian](../../../functions/geometries/geometry/DeterminantOfJacobian.md)
- [DeterminantOfJacobian](../../../functions/geometries/nurbs_surface_on_volume_geometry/NurbsSurfaceOnVolumeGeometry/DeterminantOfJacobian.md)
- [DeterminantOfJacobian](../../../functions/geometries/quadrature_point_surface_in_volume_geometry/QuadraturePointSurfaceInVolumeGeometry/DeterminantOfJacobian.md)
- [DeterminantOfJacobian](../../../functions/geometries/surface_in_nurbs_volume_geometry/SurfaceInNurbsVolumeGeometry/DeterminantOfJacobian.md)
- [Execute](../../../functions/processes/calculate_nodal_area_process/CalculateNodalAreaProcess<THistorical>/Execute.md)
- [KRATOS_TEST_CASE_IN_SUITE](../../../functions/tests/cpp_tests/geometries/test_non_rectangular_jacobian/KRATOS_TEST_CASE_IN_SUITE.md)
- [KRATOS_TEST_CASE_IN_SUITE](../../../functions/tests/cpp_tests/utilities/test_math_utils/KRATOS_TEST_CASE_IN_SUITE.md)