---
type: C++ Method
title: Weight
resource: integration/integration_point.h#L200-L203
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/geometries/brep_curve_on_surface/BrepCurveOnSurface/Length
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/geometries/line_2d_4/Line2D4/Length
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/geometries/line_2d_5/Line2D5/Length
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/geometries/nurbs_curve_geometry/NurbsCurveGeometry/Length
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/geometries/nurbs_curve_on_surface_geometry/NurbsCurveOnSurfaceGeometry/Length
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/geometries/nurbs_surface_on_volume_geometry/NurbsSurfaceOnVolumeGeometry/Area
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/geometries/point_on_geometry/PointOnGeometry/CreateIntegrationPoints
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/geometries/surface_in_nurbs_volume_geometry/SurfaceInNurbsVolumeGeometry/Area
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/geometries/triangle_2d_10/Triangle2D10/Area
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/geometries/triangle_2d_15/Triangle2D15/Area
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/integration/integration_point_utilities/IntegrationPointUtilities/CreateIntegrationPoints1DGrid
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/integration/integration_point_utilities/IntegrationPointUtilities/IntegrationPoints1D
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/integration/integration_point_utilities/IntegrationPointUtilities/IntegrationPoints2D
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/integration/integration_point_utilities/IntegrationPointUtilities/IntegrationPoints3D
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/integration/integration_point_utilities/IntegrationPointUtilities/IntegrationPointsTriangle2D
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/integration/quadrature/IntegrationPoints
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/modeler/internals/skin_intersection/SkinIntersection/AddTriangleIntegrationPointData
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/modified_shape_functions/modified_shape_functions/ModifiedShapeFunctions/ComputeValuesOnOneSide
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/modified_shape_functions/modified_shape_functions/ModifiedShapeFunctions/ComputeFaceValuesOnOneSide
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/calculate_nodal_area_process/CalculateNodalAreaProcess<THistorical>/Execute
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/compute_nodal_gradient_process/ComputeNodalGradientProcess<TOutputHistorical>/ComputeElementalContributionsAndVolume
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/compute_nodal_normal_divergence_process/ComputeNodalNormalDivergenceProcess<THistorical>/Execute
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/integration_values_extrapolation_to_nodes_process/IntegrationValuesExtrapolationToNodesProcess/ExecuteFinalizeSolutionStep
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/integration_values_extrapolation_to_nodes_process/IntegrationValuesExtrapolationToNodesProcess/InitializeMaps
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/simple_mortar_mapper_process/SimpleMortarMapperProcess<TDim, TNumNodes, TVarType, TNumNodesMaster>/AssemblyMortarOperators
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/simple_mortar_mapper_process/SimpleMortarMapperProcess<TDim, TNumNodes, TVarType, TNumNodesMaster>/CalculateAe
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tensor_adaptors/geometries_tensor_adaptor/GeometriesTensorAdaptor/CollectIntegrationWeights
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/geometries/test_brep_surface_trimmed/KRATOS_TEST_CASE_IN_SUITE
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/geometries/test_coupling_geometry_quadrature/KRATOS_TEST_CASE_IN_SUITE
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/geometries/test_geometry/CalculateAreaByIntegration
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/geometries/test_nurbs_curve/KRATOS_TEST_CASE_IN_SUITE
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/geometries/test_nurbs_curve_on_surface_geometry/KRATOS_TEST_CASE_IN_SUITE
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/geometries/test_nurbs_surface/KRATOS_TEST_CASE_IN_SUITE
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/geometries/test_nurbs_volume/KRATOS_TEST_CASE_IN_SUITE
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/geometries/test_prism_3d_15/ShapeFunctionSanityCheck
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/geometries/test_quadrature_point_geometry/KRATOS_TEST_CASE_IN_SUITE
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/geometries/test_surface_in_nurbs_volume/KRATOS_TEST_CASE_IN_SUITE
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/integration/test_integration_quadratures/KRATOS_TEST_CASE_IN_SUITE
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/sources/test_serializer_geometry_data/KRATOS_TEST_CASE_IN_SUITE
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/utilities/test_line_geometrical_sensitivity/IntegrationPointWeight
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/utilities/test_line_geometrical_sensitivity/CheckIntegrationPointWeightSensitivity
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_utilities/test_laplacian_element/TestLaplacianElement/CalculateLocalSystem
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/exact_mortar_segmentation_utility/ExactMortarIntegrationUtility<TDim, TNumNodes, TBelong, TNumNodesMaster>/GetExactIntegration
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/exact_mortar_segmentation_utility/ExactMortarIntegrationUtility<TDim, TNumNodes, TBelong, TNumNodesMaster>/TestGetExactIntegration
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/geometry_tester/GeometryTesterUtility/VerifyAreaByIntegration
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/integration_utilities/ComputeDomainSize
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/integration_utilities/ComputeArea2DGeometry
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/integration_utilities/ComputeArea2DGeometryDerivative
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/integration_utilities/ComputeVolume3DGeometry
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/integration_utilities/ComputeVolume3DGeometryDerivative
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`TWeightType Weight() const`

# Called by

- [Length](../../../../functions/geometries/brep_curve_on_surface/BrepCurveOnSurface/Length.md)
- [Length](../../../../functions/geometries/line_2d_4/Line2D4/Length.md)
- [Length](../../../../functions/geometries/line_2d_5/Line2D5/Length.md)
- [Length](../../../../functions/geometries/nurbs_curve_geometry/NurbsCurveGeometry/Length.md)
- [Length](../../../../functions/geometries/nurbs_curve_on_surface_geometry/NurbsCurveOnSurfaceGeometry/Length.md)
- [Area](../../../../functions/geometries/nurbs_surface_on_volume_geometry/NurbsSurfaceOnVolumeGeometry/Area.md)
- [CreateIntegrationPoints](../../../../functions/geometries/point_on_geometry/PointOnGeometry/CreateIntegrationPoints.md)
- [Area](../../../../functions/geometries/surface_in_nurbs_volume_geometry/SurfaceInNurbsVolumeGeometry/Area.md)
- [Area](../../../../functions/geometries/triangle_2d_10/Triangle2D10/Area.md)
- [Area](../../../../functions/geometries/triangle_2d_15/Triangle2D15/Area.md)
- [CreateIntegrationPoints1DGrid](../../../../functions/integration/integration_point_utilities/IntegrationPointUtilities/CreateIntegrationPoints1DGrid.md)
- [IntegrationPoints1D](../../../../functions/integration/integration_point_utilities/IntegrationPointUtilities/IntegrationPoints1D.md)
- [IntegrationPoints2D](../../../../functions/integration/integration_point_utilities/IntegrationPointUtilities/IntegrationPoints2D.md)
- [IntegrationPoints3D](../../../../functions/integration/integration_point_utilities/IntegrationPointUtilities/IntegrationPoints3D.md)
- [IntegrationPointsTriangle2D](../../../../functions/integration/integration_point_utilities/IntegrationPointUtilities/IntegrationPointsTriangle2D.md)
- [IntegrationPoints](../../../../functions/integration/quadrature/IntegrationPoints.md)
- [AddTriangleIntegrationPointData](../../../../functions/modeler/internals/skin_intersection/SkinIntersection/AddTriangleIntegrationPointData.md)
- [ComputeValuesOnOneSide](../../../../functions/modified_shape_functions/modified_shape_functions/ModifiedShapeFunctions/ComputeValuesOnOneSide.md)
- [ComputeFaceValuesOnOneSide](../../../../functions/modified_shape_functions/modified_shape_functions/ModifiedShapeFunctions/ComputeFaceValuesOnOneSide.md)
- [Execute](../../../../functions/processes/calculate_nodal_area_process/CalculateNodalAreaProcess<THistorical>/Execute.md)
- [ComputeElementalContributionsAndVolume](../../../../functions/processes/compute_nodal_gradient_process/ComputeNodalGradientProcess<TOutputHistorical>/ComputeElementalContributionsAndVolume.md)
- [Execute](../../../../functions/processes/compute_nodal_normal_divergence_process/ComputeNodalNormalDivergenceProcess<THistorical>/Execute.md)
- [ExecuteFinalizeSolutionStep](../../../../functions/processes/integration_values_extrapolation_to_nodes_process/IntegrationValuesExtrapolationToNodesProcess/ExecuteFinalizeSolutionStep.md)
- [InitializeMaps](../../../../functions/processes/integration_values_extrapolation_to_nodes_process/IntegrationValuesExtrapolationToNodesProcess/InitializeMaps.md)
- [AssemblyMortarOperators](../../../../functions/processes/simple_mortar_mapper_process/SimpleMortarMapperProcess<TDim, TNumNodes, TVarType, TNumNodesMaster>/AssemblyMortarOperators.md)
- [CalculateAe](../../../../functions/processes/simple_mortar_mapper_process/SimpleMortarMapperProcess<TDim, TNumNodes, TVarType, TNumNodesMaster>/CalculateAe.md)
- [CollectIntegrationWeights](../../../../functions/tensor_adaptors/geometries_tensor_adaptor/GeometriesTensorAdaptor/CollectIntegrationWeights.md)
- [KRATOS_TEST_CASE_IN_SUITE](../../../../functions/tests/cpp_tests/geometries/test_brep_surface_trimmed/KRATOS_TEST_CASE_IN_SUITE.md)
- [KRATOS_TEST_CASE_IN_SUITE](../../../../functions/tests/cpp_tests/geometries/test_coupling_geometry_quadrature/KRATOS_TEST_CASE_IN_SUITE.md)
- [CalculateAreaByIntegration](../../../../functions/tests/cpp_tests/geometries/test_geometry/CalculateAreaByIntegration.md)
- [KRATOS_TEST_CASE_IN_SUITE](../../../../functions/tests/cpp_tests/geometries/test_nurbs_curve/KRATOS_TEST_CASE_IN_SUITE.md)
- [KRATOS_TEST_CASE_IN_SUITE](../../../../functions/tests/cpp_tests/geometries/test_nurbs_curve_on_surface_geometry/KRATOS_TEST_CASE_IN_SUITE.md)
- [KRATOS_TEST_CASE_IN_SUITE](../../../../functions/tests/cpp_tests/geometries/test_nurbs_surface/KRATOS_TEST_CASE_IN_SUITE.md)
- [KRATOS_TEST_CASE_IN_SUITE](../../../../functions/tests/cpp_tests/geometries/test_nurbs_volume/KRATOS_TEST_CASE_IN_SUITE.md)
- [ShapeFunctionSanityCheck](../../../../functions/tests/cpp_tests/geometries/test_prism_3d_15/ShapeFunctionSanityCheck.md)
- [KRATOS_TEST_CASE_IN_SUITE](../../../../functions/tests/cpp_tests/geometries/test_quadrature_point_geometry/KRATOS_TEST_CASE_IN_SUITE.md)
- [KRATOS_TEST_CASE_IN_SUITE](../../../../functions/tests/cpp_tests/geometries/test_surface_in_nurbs_volume/KRATOS_TEST_CASE_IN_SUITE.md)
- [KRATOS_TEST_CASE_IN_SUITE](../../../../functions/tests/cpp_tests/integration/test_integration_quadratures/KRATOS_TEST_CASE_IN_SUITE.md)
- [KRATOS_TEST_CASE_IN_SUITE](../../../../functions/tests/cpp_tests/sources/test_serializer_geometry_data/KRATOS_TEST_CASE_IN_SUITE.md)
- [IntegrationPointWeight](../../../../functions/tests/cpp_tests/utilities/test_line_geometrical_sensitivity/IntegrationPointWeight.md)
- [CheckIntegrationPointWeightSensitivity](../../../../functions/tests/cpp_tests/utilities/test_line_geometrical_sensitivity/CheckIntegrationPointWeightSensitivity.md)
- [CalculateLocalSystem](../../../../functions/tests/test_utilities/test_laplacian_element/TestLaplacianElement/CalculateLocalSystem.md)
- [GetExactIntegration](../../../../functions/utilities/exact_mortar_segmentation_utility/ExactMortarIntegrationUtility<TDim, TNumNodes, TBelong, TNumNodesMaster>/GetExactIntegration.md)
- [TestGetExactIntegration](../../../../functions/utilities/exact_mortar_segmentation_utility/ExactMortarIntegrationUtility<TDim, TNumNodes, TBelong, TNumNodesMaster>/TestGetExactIntegration.md)
- [VerifyAreaByIntegration](../../../../functions/utilities/geometry_tester/GeometryTesterUtility/VerifyAreaByIntegration.md)
- [ComputeDomainSize](../../../../functions/utilities/integration_utilities/ComputeDomainSize.md)
- [ComputeArea2DGeometry](../../../../functions/utilities/integration_utilities/ComputeArea2DGeometry.md)
- [ComputeArea2DGeometryDerivative](../../../../functions/utilities/integration_utilities/ComputeArea2DGeometryDerivative.md)
- [ComputeVolume3DGeometry](../../../../functions/utilities/integration_utilities/ComputeVolume3DGeometry.md)
- [ComputeVolume3DGeometryDerivative](../../../../functions/utilities/integration_utilities/ComputeVolume3DGeometryDerivative.md)