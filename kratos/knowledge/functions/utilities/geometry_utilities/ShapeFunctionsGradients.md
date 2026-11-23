---
type: C++ Function
title: ShapeFunctionsGradients
resource: utilities/geometry_utilities.h#L670-L680
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/includes/amatrix_interface/prod
    resolved_by: tree-sitter
    confidence: exact
  called_by:
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

`static void ShapeFunctionsGradients( TMatrix1 const& rDN_De, TMatrix2 const& rInvJ, TMatrix3& rDN_DX )`

# Calls

- [prod](../../../functions/includes/amatrix_interface/prod.md)

# Called by

- [ComputeElementalContributionsAndVolume](../../../functions/processes/compute_nodal_gradient_process/ComputeNodalGradientProcess<TOutputHistorical>/ComputeElementalContributionsAndVolume.md)
- [Execute](../../../functions/processes/compute_nodal_normal_divergence_process/ComputeNodalNormalDivergenceProcess<THistorical>/Execute.md)
- [KRATOS_TEST_CASE_IN_SUITE](../../../functions/tests/cpp_tests/utilities/test_geometry_utils/KRATOS_TEST_CASE_IN_SUITE.md)
- [CalculateLocalSystem](../../../functions/tests/test_utilities/test_laplacian_element/TestLaplacianElement/CalculateLocalSystem.md)