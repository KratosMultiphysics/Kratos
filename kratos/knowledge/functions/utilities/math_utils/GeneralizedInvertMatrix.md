---
type: C++ Function
title: GeneralizedInvertMatrix
resource: utilities/math_utils.h#L290-L321
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/includes/amatrix_interface/prod
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/geometries/geometry/ShapeFunctionsIntegrationPointsGradients
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/compute_nodal_gradient_process/ComputeNodalGradientProcess<TOutputHistorical>/ComputeElementalContributionsAndVolume
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/compute_nodal_normal_divergence_process/ComputeNodalNormalDivergenceProcess<THistorical>/Execute
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/utilities/test_math_utils/KRATOS_TEST_CASE_IN_SUITE
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`static void GeneralizedInvertMatrix( const TMatrix1& rInputMatrix, TMatrix2& rInvertedMatrix, double& rInputMatrixDet, const double Tolerance = ZeroTolerance )`

# Calls

- [prod](../../../functions/includes/amatrix_interface/prod.md)

# Called by

- [ShapeFunctionsIntegrationPointsGradients](../../../functions/geometries/geometry/ShapeFunctionsIntegrationPointsGradients.md)
- [ComputeElementalContributionsAndVolume](../../../functions/processes/compute_nodal_gradient_process/ComputeNodalGradientProcess<TOutputHistorical>/ComputeElementalContributionsAndVolume.md)
- [Execute](../../../functions/processes/compute_nodal_normal_divergence_process/ComputeNodalNormalDivergenceProcess<THistorical>/Execute.md)
- [KRATOS_TEST_CASE_IN_SUITE](../../../functions/tests/cpp_tests/utilities/test_math_utils/KRATOS_TEST_CASE_IN_SUITE.md)