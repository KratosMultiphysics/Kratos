---
type: C++ Function
title: JacobianOnInitialConfiguration
resource: utilities/geometry_utilities.h#L708-L720
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/geometries/geometry/PointsNumber
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/processes/calculate_nodal_area_process/CalculateNodalAreaProcess<THistorical>/Execute
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/compute_nodal_gradient_process/ComputeNodalGradientProcess<TOutputHistorical>/ComputeElementalContributionsAndVolume
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/compute_nodal_normal_divergence_process/ComputeNodalNormalDivergenceProcess<THistorical>/Execute
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/utilities/test_geometry_utils/KRATOS_TEST_CASE_IN_SUITE
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_utilities/test_laplacian_element/TestLaplacianElement/CalculateLocalSystem
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`static void JacobianOnInitialConfiguration( GeometryType const& rGeom, GeometryType::CoordinatesArrayType const& rCoords, Matrix& rJ0 )`

# Calls

- [PointsNumber](../../../functions/geometries/geometry/PointsNumber.md)

# Called by

- [Execute](../../../functions/processes/calculate_nodal_area_process/CalculateNodalAreaProcess<THistorical>/Execute.md)
- [ComputeElementalContributionsAndVolume](../../../functions/processes/compute_nodal_gradient_process/ComputeNodalGradientProcess<TOutputHistorical>/ComputeElementalContributionsAndVolume.md)
- [Execute](../../../functions/processes/compute_nodal_normal_divergence_process/ComputeNodalNormalDivergenceProcess<THistorical>/Execute.md)
- [KRATOS_TEST_CASE_IN_SUITE](../../../functions/tests/cpp_tests/utilities/test_geometry_utils/KRATOS_TEST_CASE_IN_SUITE.md)
- [CalculateLocalSystem](../../../functions/tests/test_utilities/test_laplacian_element/TestLaplacianElement/CalculateLocalSystem.md)